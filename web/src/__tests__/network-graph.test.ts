import { describe, expect, it } from "vitest";

import { centeredOffsets, edgePath } from "../lib/network-geometry";
import type { GraphEdge, GraphNode } from "../types/cytocarto";

describe("network edge geometry", () => {
  it("stops directed edges before the target gene node", () => {
    const source: GraphNode & { x: number; y: number } = { id: "analyte:IL-6", label: "IL-6", kind: "analyte", x: 145, y: 100 };
    const target: GraphNode & { x: number; y: number } = { id: "gene:STAT3", label: "STAT3", kind: "gene", x: 335, y: 100 };
    const edge: GraphEdge = { id: "test", source: source.id, target: target.id, kind: "analyte_gene", weight: 1 };

    expect(edgePath(edge, source, target)).toBe("M 164 100 C 232 100, 232 100, 300 100");
  });

  it("centers shared edge attachment points with close, even spacing", () => {
    expect(centeredOffsets(1)).toEqual([0]);
    expect(centeredOffsets(2)).toEqual([-4, 4]);
    expect(centeredOffsets(4)).toEqual([-12, -4, 4, 12]);
  });

  it("keeps angled incoming edges centered on the target node's left side", () => {
    const source: GraphNode & { x: number; y: number } = { id: "analyte:IL-10", label: "IL-10", kind: "analyte", x: 145, y: 200 };
    const target: GraphNode & { x: number; y: number } = { id: "gene:FOS", label: "FOS", kind: "gene", x: 335, y: 100 };
    const edge: GraphEdge = { id: "angled", source: source.id, target: target.id, kind: "analyte_gene", weight: 1 };
    const path = edgePath(edge, source, target, { targetOffset: 4 });

    expect(path.endsWith(`${335 - Math.sqrt(1209)} 104`)).toBe(true);
  });
});
