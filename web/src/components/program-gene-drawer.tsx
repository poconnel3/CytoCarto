"use client";

import { useEffect, useMemo, useRef, useState } from "react";
import { Download, X } from "lucide-react";
import { useVirtualizer } from "@tanstack/react-virtual";

import { getProgramGenes } from "@/lib/api";
import type { CytoCartoResult } from "@/types/cytocarto";

type ProgramGeneDrawerProps = {
  program: string | null;
  result: CytoCartoResult | null;
  onClose: () => void;
  variant?: "drawer" | "panel";
};

export function ProgramGeneDrawer({ program, result, onClose, variant = "drawer" }: ProgramGeneDrawerProps) {
  const scrollRef = useRef<HTMLDivElement>(null);
  const [query, setQuery] = useState("");
  const [members, setMembers] = useState<Array<{ gene: string; program_weight: number; direction: string }>>([]);
  const [loading, setLoading] = useState(false);
  const geneScores = useMemo(() => new Map((result?.tables["gene_subnetwork.tsv"] ?? []).map((row) => [String(row.gene), Number(row.score ?? 0)])), [result]);
  const programRow = useMemo(() => (result?.tables["program_scores.tsv"] ?? []).find((row) => String(row.program) === program), [program, result]);
  const rows = useMemo(() => program && result ? members : [], [members, program, result]);
  useEffect(() => {
    if (!program || !result) return;
    setLoading(true);
    setQuery("");
    getProgramGenes(program).then(setMembers).catch(() => setMembers([])).finally(() => setLoading(false));
  }, [program, result]);
  const visibleRows = useMemo(() => rows.filter((row) => row.gene.toLowerCase().includes(query.trim().toLowerCase())), [query, rows]);
  // TanStack Virtual owns imperative measurements and intentionally returns methods.
  // eslint-disable-next-line react-hooks/incompatible-library
  const virtualizer = useVirtualizer({ count: visibleRows.length, getScrollElement: () => scrollRef.current, estimateSize: () => 52, overscan: 10 });

  if (!program || !result) return null;
  const programLabel = String(programRow?.final_annotation ?? `P${programRow?.program_id ?? ""}`)
    .replace(/na(?:ï|\uFFFD)ve/gi, (value) => value[0] === "N" ? "Naive" : "naive")
    .replace(/^P\d+[- ]?/, "")
    .replace("reprogramed", "reprogrammed");
  const exportGenes = () => {
    const header = "gene\tdirection\tprogram_weight\tpatient_gene_score";
    const lines = rows.map((row) => `${row.gene}\t${row.direction}\t${row.program_weight}\t${geneScores.get(row.gene) ?? 0}`);
    const url = URL.createObjectURL(new Blob([[header, ...lines].join("\n")], { type: "text/tab-separated-values" }));
    const anchor = document.createElement("a");
    anchor.href = url;
    anchor.download = `cytocarto-${String(programRow?.program_id ?? "program")}-genes.tsv`;
    anchor.click();
    URL.revokeObjectURL(url);
  };

  return (
    <aside className={`gene-drawer ${variant === "panel" ? "gene-panel" : ""}`} aria-label={`Genes associated with ${program}`}>
      <div className="drawer-header">
        <div>
          <h2><span className="program-code">P{String(programRow?.program_id ?? "")}</span> <i>—</i> {programLabel}</h2>
          <dl className="program-summary">
            <div><dt>Score</dt><dd>{Number(programRow?.score ?? 0).toFixed(3)}</dd></div>
            <div><dt>Annotation</dt><dd>{String(programRow?.biological_response_or_function ?? programRow?.immune_cell_type ?? "DSPIN program")}</dd></div>
          </dl>
        </div>
        {variant === "drawer" ? <button className="icon-button" onClick={onClose} aria-label="Close gene list"><X size={20} /></button> : null}
      </div>
      <div className="gene-table-heading"><strong>Associated genes</strong><span>({loading ? "loading" : rows.length.toLocaleString()})</span></div>
      <input className="gene-search" value={query} onChange={(event) => setQuery(event.target.value)} placeholder="Search genes" aria-label="Search associated genes" />
      <div className="gene-column-head"><span>Gene</span><span>Direction</span><span>Program weight</span><span>Patient score</span></div>
      <div ref={scrollRef} className="gene-list">
        <div style={{ height: `${virtualizer.getTotalSize()}px`, position: "relative" }}>
          {virtualizer.getVirtualItems().map((item) => {
            const row = visibleRows[item.index];
            const patientScore = geneScores.get(row.gene) ?? 0;
            return <div className="gene-row" key={row.gene} style={{ transform: `translateY(${item.start}px)` }}>
              <strong>{row.gene}</strong>
              <span className={row.direction === "positive" ? "positive" : "negative"}>{row.direction === "positive" ? "↑" : "↓"}</span>
              <span>{row.program_weight >= 0 ? "+" : ""}{row.program_weight.toFixed(3)}</span>
              <span className={patientScore >= 0 ? "positive" : "negative"}>{patientScore.toFixed(2)}</span>
            </div>;
          })}
        </div>
      </div>
      <button className="button button-secondary export-genes" onClick={exportGenes}><Download size={17} /> Export genes</button>
    </aside>
  );
}
