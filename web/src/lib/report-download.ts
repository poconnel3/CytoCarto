import type { CytoCartoResult, TableRow } from "@/types/cytocarto";

function tsvValue(value: unknown): string {
  if (value === null || value === undefined) return "";
  const text = typeof value === "object" ? JSON.stringify(value) : String(value);
  return /[\t\r\n"]/.test(text) ? `"${text.replaceAll('"', '""')}"` : text;
}

export function tableRowsToTsv(rows: TableRow[]): string {
  if (!rows.length) return "";
  const columns = Array.from(new Set(rows.flatMap((row) => Object.keys(row))));
  const lines = [columns.map(tsvValue).join("\t")];
  for (const row of rows) lines.push(columns.map((column) => tsvValue(row[column])).join("\t"));
  return `${lines.join("\n")}\n`;
}

async function browserResultBundle(result: CytoCartoResult): Promise<Blob> {
  const { default: JSZip } = await import("jszip");
  const archive = new JSZip();
  archive.file("patient_report.json", `${JSON.stringify(result.report, null, 2)}\n`);
  archive.file("cytocarto_result.json", `${JSON.stringify(result, null, 2)}\n`);
  for (const [filename, rows] of Object.entries(result.tables)) {
    archive.file(filename, tableRowsToTsv(rows));
  }
  archive.file(
    "README.txt",
    [
      "CytoCarto report export",
      "",
      "This ZIP was generated from the structured analysis result already returned to this browser.",
      "Very large evidence tables may contain only the rows returned by the analysis API.",
      "patient_report.json contains the complete report object used by the interface.",
      "",
    ].join("\n"),
  );
  return archive.generateAsync({ type: "blob", compression: "DEFLATE" });
}

function saveBlob(blob: Blob, filename: string): void {
  const url = URL.createObjectURL(blob);
  const anchor = document.createElement("a");
  anchor.href = url;
  anchor.download = filename;
  anchor.style.display = "none";
  document.body.appendChild(anchor);
  anchor.click();
  window.setTimeout(() => {
    anchor.remove();
    URL.revokeObjectURL(url);
  }, 1_000);
}

export async function downloadReport(result: CytoCartoResult, bundleUrl: string): Promise<"server" | "browser"> {
  try {
    const response = await fetch(bundleUrl);
    if (!response.ok) throw new Error(`Bundle request failed (${response.status}).`);
    const blob = await response.blob();
    if (!blob.size) throw new Error("Bundle response was empty.");
    saveBlob(blob, "cytocarto-report.zip");
    return "server";
  } catch {
    saveBlob(await browserResultBundle(result), "cytocarto-report.zip");
    return "browser";
  }
}
