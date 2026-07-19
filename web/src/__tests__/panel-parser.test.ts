import { describe, expect, it } from "vitest";

import { normalizedAnalytesForApi } from "../lib/api";
import { derivedAnalyteFlag, parseDecimalDraft } from "../lib/analyte-entry";
import { EXAMPLE_PANEL_TEXT, parsePanelText } from "../lib/panel-parser";

describe("clinical panel parser", () => {
  it("parses wrapped ARUP names, flags, censored values, and laboratory ranges", () => {
    const rows = parsePanelText(EXAMPLE_PANEL_TEXT);
    expect(rows).toHaveLength(13);
    expect(rows.find((row) => row.raw_name === "Tumor Necrosis Factor - alpha")).toMatchObject({ value: 5.4, reference_high: 7.2 });
    expect(rows.find((row) => row.raw_name === "Interleukin 2 Receptor, Soluble")).toMatchObject({ value: 4086.1, flag: "H", reference_low: 175.3, reference_high: 858.2 });
    expect(rows.find((row) => row.raw_name === "Interleukin 2")).toMatchObject({ value: 2.1, qualifier: "<", reference_high: 2.1 });
    expect(rows.find((row) => row.raw_name === "Interleukin 6")).toMatchObject({ value: 199.8, flag: "H", reference_high: 2 });
  });

  it("parses CSV and preserves a report-specific upper reference", () => {
    const rows = parsePanelText("Component,Result,Units,Ref Range\nIL-8,6.9,pg/mL,<=3.0");
    expect(rows).toEqual([expect.objectContaining({ raw_name: "IL-8", value: 6.9, reference_high: 3 })]);
  });

  it("accepts basic RTF text exports", () => {
    const rows = parsePanelText("{\\rtf1\\ansi Interleukin 5 6.8 H pg/mL <=2.1\\par}");
    expect(rows).toEqual([expect.objectContaining({ raw_name: "Interleukin 5", value: 6.8, flag: "H", reference_high: 2.1 })]);
  });

  it("keeps parser source text in-browser instead of sending it to the API", () => {
    const [row] = parsePanelText("Interleukin 6 10 H pg/mL <=2.0");
    expect(row.source).toBeTruthy();
    expect(normalizedAnalytesForApi([row])).toEqual([
      expect.not.objectContaining({ source: expect.anything() }),
    ]);
  });

  it("keeps incomplete decimal drafts editable and commits completed decimals", () => {
    expect(parseDecimalDraft(".")).toBeNull();
    expect(parseDecimalDraft("0.")).toBe(0);
    expect(parseDecimalDraft(".5")).toBe(0.5);
  });

  it("recomputes high and low flags from edited values and ranges", () => {
    const base = { raw_name: "IL-6", value: 5, units: "pg/mL" as const, flag: "H" as const, reference_low: 2, reference_high: 10 };
    expect(derivedAnalyteFlag(base)).toBeUndefined();
    expect(derivedAnalyteFlag({ ...base, value: 11 })).toBe("H");
    expect(derivedAnalyteFlag({ ...base, value: 1 })).toBe("L");
  });
});
