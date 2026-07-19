// @vitest-environment jsdom

import { afterEach, describe, expect, it, vi } from "vitest";

import { downloadReport, tableRowsToTsv } from "../lib/report-download";
import type { CytoCartoResult } from "../types/cytocarto";

afterEach(() => {
  vi.restoreAllMocks();
  vi.unstubAllGlobals();
});

describe("report download", () => {
  it("serializes returned report tables as valid tab-separated text", () => {
    expect(tableRowsToTsv([
      { rank: 1, disease_name: "Example disease", matched_genes: "STAT1;STAT3" },
      { rank: 2, disease_name: "Term with\ttab", matched_genes: null },
    ])).toBe([
      "rank\tdisease_name\tmatched_genes",
      "1\tExample disease\tSTAT1;STAT3",
      "2\t\"Term with\ttab\"\t",
      "",
    ].join("\n"));
  });

  it("downloads a browser-generated ZIP when the server bundle is unavailable", async () => {
    vi.stubGlobal("fetch", vi.fn().mockResolvedValue({ ok: false, status: 404 }));
    const createObjectUrl = vi.fn().mockReturnValue("blob:cytocarto-report");
    vi.stubGlobal("URL", { ...URL, createObjectURL: createObjectUrl, revokeObjectURL: vi.fn() });
    const click = vi.spyOn(HTMLAnchorElement.prototype, "click").mockImplementation(() => undefined);
    const result = {
      report: {
        patient: {},
        analytes: [],
        warnings: [],
        notes: [],
        cell_type_enrichment_mode: "exact_per_cell",
      },
      tables: { "human_disease_hypotheses.tsv": [{ rank: 1, disease_name: "Example disease" }] },
      download_url: "/v1/bundles/expired",
      graph: { nodes: [], edges: [] },
    } satisfies CytoCartoResult;

    await expect(downloadReport(result, "http://127.0.0.1:8000/v1/bundles/expired")).resolves.toBe("browser");
    expect(createObjectUrl).toHaveBeenCalledWith(expect.any(Blob));
    expect(click).toHaveBeenCalledOnce();
  });
});
