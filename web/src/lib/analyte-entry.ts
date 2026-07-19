import type { ParsedAnalyte } from "@/types/cytocarto";

export function parseDecimalDraft(value: string): number | null {
  const trimmed = value.trim();
  if (!/^(?:\d+(?:\.\d*)?|\.\d+)$/.test(trimmed)) return null;
  const parsed = Number(trimmed);
  return Number.isFinite(parsed) && parsed >= 0 ? parsed : null;
}

export function derivedAnalyteFlag(analyte: ParsedAnalyte): "H" | "L" | undefined {
  if (analyte.reference_high !== undefined && analyte.value > analyte.reference_high) return "H";
  if (analyte.reference_low !== undefined && analyte.value < analyte.reference_low) return "L";
  if (analyte.reference_low !== undefined || analyte.reference_high !== undefined) return undefined;
  return analyte.flag;
}
