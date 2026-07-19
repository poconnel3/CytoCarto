"use client";

import { useMemo } from "react";
import { Dna } from "lucide-react";

import { ImmuneCellIcon } from "@/components/immune-cell-icon";
import { centeredOffsets, edgePath, NETWORK_ARROW_SIZE } from "@/lib/network-geometry";
import type { EdgeRoute } from "@/lib/network-geometry";
import type { CytoCartoResult, GraphEdge, GraphNode } from "@/types/cytocarto";

type PositionedNode = GraphNode & { x: number; y: number };

const LAYER_X: Record<GraphNode["kind"], number> = {
  analyte: 145,
  gene: 335,
  program: 565,
  cell: 770,
};

function distribute(nodes: GraphNode[], kind: GraphNode["kind"]): PositionedNode[] {
  const top = 76;
  const bottom = 670;
  return nodes.map((node, index) => ({
    ...node,
    x: LAYER_X[kind],
    y: nodes.length === 1 ? (top + bottom) / 2 : top + index * ((bottom - top) / (nodes.length - 1)),
  }));
}

function ProgramGlyph() {
  return <g aria-hidden="true"><circle r="16" className="program-cell" /><circle r="6" className="program-nucleus" /><circle cx="-11" cy="-9" r="1.7" /><circle cx="12" cy="7" r="1.7" /><circle cx="9" cy="-11" r="1.3" /></g>;
}

function analyteLabel(label: string): string {
  return label.toLowerCase().includes("soluble il-2") ? "sIL-2R alpha" : label;
}

export function NetworkGraph({ graph, selectedProgram, onProgramSelect }: { graph: CytoCartoResult["graph"]; selectedProgram: string | null; onProgramSelect: (program: string) => void }) {
  const { nodes, edges, routes } = useMemo(() => {
    const byKind = {
      analyte: graph.nodes.filter((node) => node.kind === "analyte").slice(0, 6),
      gene: graph.nodes.filter((node) => node.kind === "gene").slice(0, 12),
      program: graph.nodes.filter((node) => node.kind === "program").slice(0, 5),
      cell: graph.nodes.filter((node) => node.kind === "cell"),
    };
    const positioned = [
      ...distribute(byKind.analyte, "analyte"),
      ...distribute(byKind.gene, "gene"),
      ...distribute(byKind.program, "program"),
      ...distribute(byKind.cell, "cell"),
    ];
    const visible = new Set(positioned.map((node) => node.id));
    const seenEdges = new Set<string>();
    const seenProgramPairs = new Set<string>();
    const visibleEdges = graph.edges
      .filter((edge) => visible.has(edge.source) && visible.has(edge.target))
      .filter((edge) => {
        if (seenEdges.has(edge.id)) return false;
        seenEdges.add(edge.id);
        if (edge.kind !== "program_program") return true;
        const pair = [edge.source, edge.target].sort().join("::");
        if (seenProgramPairs.has(pair)) return false;
        seenProgramPairs.add(pair);
        return true;
      });
    const positionedById = new Map(positioned.map((node) => [node.id, node]));
    const routesByEdge = new Map<string, EdgeRoute>();
    const incoming = new Map<string, GraphEdge[]>();
    const outgoing = new Map<string, GraphEdge[]>();
    for (const edge of visibleEdges) {
      if (edge.kind !== "analyte_gene" && edge.kind !== "gene_program") continue;
      incoming.set(edge.target, [...(incoming.get(edge.target) ?? []), edge]);
      outgoing.set(edge.source, [...(outgoing.get(edge.source) ?? []), edge]);
    }
    for (const group of incoming.values()) {
      group.sort((left, right) => (positionedById.get(left.source)?.y ?? 0) - (positionedById.get(right.source)?.y ?? 0));
      centeredOffsets(group.length).forEach((offset, index) => {
        routesByEdge.set(group[index].id, { ...(routesByEdge.get(group[index].id) ?? {}), targetOffset: offset });
      });
    }
    for (const group of outgoing.values()) {
      group.sort((left, right) => (positionedById.get(left.target)?.y ?? 0) - (positionedById.get(right.target)?.y ?? 0));
      centeredOffsets(group.length).forEach((offset, index) => {
        routesByEdge.set(group[index].id, { ...(routesByEdge.get(group[index].id) ?? {}), sourceOffset: offset });
      });
    }
    return {
      nodes: positioned,
      edges: visibleEdges,
      routes: routesByEdge,
    };
  }, [graph]);
  const byId = new Map(nodes.map((node) => [node.id, node]));
  const maxContextWeight = Math.max(
    ...edges.filter((edge) => edge.kind === "cell_context").map((edge) => Math.abs(edge.weight)),
    0,
  );

  return <div className="network-diagram" aria-label="Signed cytokine to DSPIN network">
    <svg viewBox="0 0 960 740" preserveAspectRatio="xMidYMid meet" role="img" aria-labelledby="network-title network-description">
      <title id="network-title">Cytokine response and DSPIN program network</title>
      <desc id="network-description">Measured analytes connect to response genes, genes connect to DSPIN programs, programs connect directionally, and dashed links show cell context.</desc>
      <defs>
        <marker id="arrow-positive" markerUnits="userSpaceOnUse" markerWidth={NETWORK_ARROW_SIZE} markerHeight={NETWORK_ARROW_SIZE} viewBox="0 0 11 11" refX="0.5" refY="5.5" orient="auto"><path d="M0.5 1 10.5 5.5 0.5 10Z" /></marker>
        <marker id="arrow-negative" markerUnits="userSpaceOnUse" markerWidth={NETWORK_ARROW_SIZE} markerHeight={NETWORK_ARROW_SIZE} viewBox="0 0 11 11" refX="0.5" refY="5.5" orient="auto"><path d="M0.5 1 10.5 5.5 0.5 10Z" /></marker>
      </defs>
      <g className="network-edges">
        {edges.map((edge) => {
          const source = byId.get(edge.source);
          const target = byId.get(edge.target);
          if (!source || !target) return null;
          const negative = edge.weight < 0;
          const contextWidth = maxContextWeight > 0 ? 2 + 4 * (Math.abs(edge.weight) / maxContextWeight) : 2;
          return <path
            key={edge.id}
            d={edgePath(edge, source, target, routes.get(edge.id))}
            className={`${edge.kind === "cell_context" ? "context-edge" : negative ? "negative-edge" : "positive-edge"}`}
            style={{ strokeWidth: edge.kind === "cell_context" ? contextWidth : 3.8 }}
            markerEnd={edge.kind === "cell_context" ? undefined : negative ? "url(#arrow-negative)" : "url(#arrow-positive)"}
          >{edge.kind === "cell_context" ? <title>{`Mean cosine contribution: ${edge.weight.toFixed(4)}`}</title> : null}</path>;
        })}
      </g>
      <g className="network-nodes">
        {nodes.map((node) => {
          if (node.kind === "analyte") return <g key={node.id} transform={`translate(${node.x} ${node.y})`} className="analyte-node">
            <text x="-31" y="5" textAnchor="end">{analyteLabel(node.label)}</text><circle r="19" /><path d="M0 10V-9M-7-2 0-9l7 7" />
          </g>;
          if (node.kind === "gene") return <g key={node.id} transform={`translate(${node.x} ${node.y})`} className="gene-node">
            <circle r="22" /><Dna x={-14} y={-14} width={28} height={28} strokeWidth={1.9} aria-hidden="true" /><text x="35" y="5">{node.label}</text>
          </g>;
          if (node.kind === "program") {
            const program = node.id.replace(/^program:/, "");
            const selected = selectedProgram === program;
            return <g
              key={node.id}
              transform={`translate(${node.x} ${node.y})`}
              className={`program-node${selected ? " selected" : ""}`}
              role="button"
              tabIndex={0}
              aria-label={`Open genes for ${node.label}`}
              onClick={() => onProgramSelect(program)}
              onKeyDown={(event) => { if (event.key === "Enter" || event.key === " ") onProgramSelect(program); }}
            ><circle r="31" /><ProgramGlyph /><text y="49" textAnchor="middle">{node.label}</text></g>;
          }
          return <g key={node.id} transform={`translate(${node.x} ${node.y})`} className="cell-node">
            <circle r="32" /><ImmuneCellIcon cellType={node.label} x={-23} y={-23} size={46} className="network-cell-symbol" /><text x="44" y="5">{node.label}</text>
          </g>;
        })}
      </g>
    </svg>
  </div>;
}
