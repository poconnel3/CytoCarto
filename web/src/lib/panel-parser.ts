import Papa from "papaparse";

import type { ParsedAnalyte } from "@/types/cytocarto";

const UNIT_FACTORS: Record<string, ParsedAnalyte["units"]> = {
  "pg/ml": "pg/mL",
  "ng/l": "ng/L",
  "ng/ml": "ng/mL",
};

function numberFrom(value: string | number | undefined): number | undefined {
  if (value === undefined || value === null || String(value).trim() === "") return undefined;
  const parsed = Number(String(value).replace(/,/g, "").replace(/^[<>]=?/, "").trim());
  return Number.isFinite(parsed) ? parsed : undefined;
}

function qualifierFrom(value: string | number | undefined): ParsedAnalyte["qualifier"] {
  const text = String(value ?? "").trim();
  return text.startsWith("<") ? "<" : text.startsWith(">") ? ">" : "=";
}

function parseReference(value: string | number | undefined): Pick<ParsedAnalyte, "reference_low" | "reference_high"> {
  const text = String(value ?? "").trim().replace(/,/g, "");
  if (!text) return {};
  const range = text.match(/(\d+(?:\.\d+)?)\s*-\s*(\d+(?:\.\d+)?)/);
  if (range) return { reference_low: Number(range[1]), reference_high: Number(range[2]) };
  const upper = text.match(/(?:<=|<)\s*(\d+(?:\.\d+)?)/);
  if (upper) return { reference_high: Number(upper[1]) };
  const lower = text.match(/(?:>=|>)\s*(\d+(?:\.\d+)?)/);
  if (lower) return { reference_low: Number(lower[1]) };
  const exact = numberFrom(text);
  return exact === undefined ? {} : { reference_high: exact };
}

function sanitizeName(value: string): string {
  return value.replace(/\s+/g, " ").replace(/\s+,/g, ",").trim();
}

function parseRow(row: Record<string, unknown>): ParsedAnalyte | null {
  const find = (...names: string[]) => {
    const key = Object.keys(row).find((candidate) => names.includes(candidate.trim().toLowerCase()));
    return key ? row[key] : undefined;
  };
  const rawName = find("analyte", "component", "test", "name", "cytokine");
  const rawValue = find("value", "result", "results", "concentration");
  const name = sanitizeName(String(rawName ?? ""));
  const value = numberFrom(rawValue as string);
  if (!name || value === undefined) return null;
  const unitsRaw = String(find("units", "unit") ?? "pg/mL").replace(/\s/g, "").toLowerCase();
  const units = UNIT_FACTORS[unitsRaw] ?? "pg/mL";
  const flag = String(find("flag", "abnormal") ?? "").trim().toUpperCase();
  return {
    raw_name: name,
    value,
    qualifier: qualifierFrom(rawValue as string),
    units,
    flag: flag === "H" || flag === "L" ? flag : undefined,
    ...parseReference(find("ref range", "reference", "reference range", "refvalue") as string),
  };
}

function stripRtf(input: string): string {
  return input
    .replace(/\\par[d]?/g, "\n")
    .replace(/\\'[0-9a-fA-F]{2}/g, "")
    .replace(/\\[a-z]+-?\d* ?/g, "")
    .replace(/[{}]/g, "");
}

function readTextLines(input: string): ParsedAnalyte[] {
  const lines = stripRtf(input).replace(/\r/g, "").split("\n");
  const parsed: ParsedAnalyte[] = [];
  const rowPattern = /^\s*([A-Za-z][A-Za-z0-9,()\- /]+?)\s+([<>]=?\s*)?([\d,.]+)\s*([HL])?\s*(pg\/mL|ng\/L|ng\/mL)\s*(.*)$/i;
  for (let index = 0; index < lines.length; index += 1) {
    const match = lines[index].match(rowPattern);
    if (!match) continue;
    let name = sanitizeName(match[1]);
    const continuation = lines[index + 1]?.trim();
    if (continuation && /^(alpha|beta|gamma|soluble|receptor|a|b)$/i.test(continuation)) {
      name = `${name} ${continuation}`;
      index += 1;
    }
    const value = numberFrom(match[3]);
    if (value === undefined) continue;
    const unitKey = match[5].replace(/\s/g, "").toLowerCase();
    parsed.push({
      raw_name: name,
      value,
      qualifier: match[2]?.trim().startsWith("<") ? "<" : match[2]?.trim().startsWith(">") ? ">" : "=",
      units: UNIT_FACTORS[unitKey] ?? "pg/mL",
      flag: match[4]?.toUpperCase() === "H" || match[4]?.toUpperCase() === "L" ? (match[4].toUpperCase() as "H" | "L") : undefined,
      ...parseReference(match[6]),
      source: lines[index],
    });
  }
  return parsed;
}

function parseJson(input: string): ParsedAnalyte[] | null {
  try {
    const payload = JSON.parse(input) as unknown;
    const rows = Array.isArray(payload) ? payload : typeof payload === "object" && payload !== null ? (payload as { analytes?: unknown[]; cytokines?: Record<string, unknown> }) : null;
    if (Array.isArray(rows)) return rows.map((row) => parseRow(row as Record<string, unknown>)).filter(Boolean) as ParsedAnalyte[];
    if (rows?.analytes) return rows.analytes.map((row) => parseRow(row as Record<string, unknown>)).filter(Boolean) as ParsedAnalyte[];
    if (rows?.cytokines) return Object.entries(rows.cytokines).map(([raw_name, value]) => ({ raw_name, value: numberFrom(value as string) ?? 0, qualifier: qualifierFrom(value as string), units: "pg/mL" as const }));
  } catch {
    return null;
  }
  return null;
}

export function parsePanelText(input: string): ParsedAnalyte[] {
  const jsonRows = parseJson(input);
  if (jsonRows?.length) return jsonRows;
  const firstLine = input.split(/\r?\n/).find((line) => line.trim()) ?? "";
  if (firstLine.includes(",") || firstLine.includes("\t")) {
    const result = Papa.parse<Record<string, unknown>>(input, { header: true, skipEmptyLines: true });
    const rows = result.data.map(parseRow).filter(Boolean) as ParsedAnalyte[];
    if (rows.length) return rows;
  }
  return readTextLines(input);
}

export async function parsePanelFile(file: File): Promise<ParsedAnalyte[]> {
  const lowerName = file.name.toLowerCase();
  if (lowerName.endsWith(".pdf")) {
    const pdfjs = await import("pdfjs-dist/legacy/build/pdf.mjs");
    const document = await pdfjs.getDocument({ data: new Uint8Array(await file.arrayBuffer()) }).promise;
    let text = "";
    for (let pageNumber = 1; pageNumber <= document.numPages; pageNumber += 1) {
      const page = await document.getPage(pageNumber);
      const content = await page.getTextContent();
      text += `${content.items.map((item) => ("str" in item ? item.str : "")).join(" ")}\n`;
    }
    const parsed = parsePanelText(text);
    if (!parsed.length) throw new Error("No extractable cytokine rows were found. Scanned PDFs require manual entry or pasted text.");
    return parsed;
  }
  return parsePanelText(await file.text());
}

export const EXAMPLE_PANEL_TEXT = `Tumor Necrosis Factor -     5.4           pg/mL  <=7.2
  alpha
  Interleukin 2               <2.1          pg/mL  <=2.1
  Interleukin 2 Receptor,     4086.1   H    pg/mL  175.3-858.2
    Soluble
  Interleukin 12              <1.9          pg/mL  <=1.9
  Interferon gamma            <4.2          pg/mL  <=4.2
  Interleukin 4               <2.2          pg/mL  <=2.2
  Interleukin 5               6.8      H    pg/mL  <=2.1
  Interleukin 10              186.7    H    pg/mL  <=2.8
  Interleukin 13              3.1      H    pg/mL  <=2.3
  Interleukin 17              <1.4          pg/mL  <=1.4
  Interleukin 1 beta          9.0      H    pg/mL  <=6.7
  Interleukin 6               199.8    H    pg/mL  <=2.0
  Interleukin 8               6.9      H    pg/mL  <=3.0`;
