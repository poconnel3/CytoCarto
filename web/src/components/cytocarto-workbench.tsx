"use client";

import { ChangeEvent, Fragment, useCallback, useEffect, useMemo, useState } from "react";
import type { CSSProperties, KeyboardEvent as ReactKeyboardEvent, PointerEvent as ReactPointerEvent } from "react";
import { MarkGithubIcon } from "@primer/octicons-react";
import {
  Check,
  CircleHelp,
  ClipboardPaste,
  Dna,
  Download,
  FileText,
  FlaskConical,
  LoaderCircle,
  Map as MapIcon,
  Network,
  Pencil,
  Plus,
  Search,
  TableProperties,
  Trash2,
  Upload,
  X,
} from "lucide-react";

import { getNetworkCytokines, apiPath, runAnalysis } from "@/lib/api";
import { derivedAnalyteFlag, parseDecimalDraft } from "@/lib/analyte-entry";
import { EXAMPLE_PANEL_TEXT, parsePanelFile, parsePanelText } from "@/lib/panel-parser";
import { downloadReport } from "@/lib/report-download";
import { ImmuneCellIcon } from "@/components/immune-cell-icon";
import { NetworkGraph } from "@/components/network-graph";
import { ProgramGeneDrawer } from "@/components/program-gene-drawer";
import type { CytoCartoResult, NetworkCytokine, ParsedAnalyte, TableRow } from "@/types/cytocarto";

type Tab = "Overview" | "Network" | "Programs" | "Cell types" | "Genes" | "Disease & mimics" | "Evidence";
type InputMode = "paste" | "upload" | "manual";

const tabs: Tab[] = ["Overview", "Network", "Programs", "Cell types", "Genes", "Disease & mimics", "Evidence"];

const RACE_OPTIONS = [
  ["American Indian or Alaska Native", "American Indian or Alaska Native"],
  ["Asian", "Asian"],
  ["Black or African American", "Black or African American"],
  ["Middle Eastern or North African", "Middle Eastern or North African"],
  ["Native Hawaiian or Pacific Islander", "Native Hawaiian or Pacific Islander"],
  ["White / Caucasian", "Caucasian"],
  ["Multiracial", "Multiracial"],
  ["Not reported", "Not reported"],
] as const;

const ETHNICITY_OPTIONS = ["Hispanic/Latino", "Not Hispanic/Latino", "Not reported"] as const;

const EVIDENCE_CHANNELS = [
  ["direct_score", "Direct"],
  ["dspin_score", "DSPIN"],
  ["chea_score", "ChEA"],
  ["manifold_score", "Manifold"],
  ["signed_sensitivity_score", "Signed sensitivity"],
  ["ppi_score", "PPI"],
  ["coexpression_score", "Co-expression"],
] as const;

type NumericAnalyteKey = "value" | "reference_low" | "reference_high";

function withEntryIds(rows: ParsedAnalyte[], manual = false): ParsedAnalyte[] {
  return rows.map((row) => ({ ...row, entry_id: crypto.randomUUID(), manual }));
}

function draftKey(row: ParsedAnalyte, index: number, field: NumericAnalyteKey): string {
  return `${row.entry_id ?? index}:${field}`;
}

function asRows(result: CytoCartoResult | null, name: string): TableRow[] {
  return result?.tables[name] ?? [];
}

function cleanDisplayText(value: unknown): string {
  return String(value ?? "")
    .replace(/na(?:ï|\uFFFD)ve/gi, (match) => match[0] === "N" ? "Naive" : "naive")
    .replace("reprogramed", "reprogrammed");
}

function compactProgramLabel(row: TableRow): string {
  return cleanDisplayText(row.final_annotation ?? row.program ?? "DSPIN program")
    .replace(/^P\d+[- ]?/, "")
    .trim();
}

function ResultTable({ rows, columns }: { rows: TableRow[]; columns: Array<{ key: string; label: string; format?: (value: unknown) => string }> }) {
  if (!rows.length) return <p className="empty-state">No result rows are available for this view.</p>;
  return <div className="data-table-wrap"><table className="data-table"><thead><tr>{columns.map((column) => <th key={column.key}>{column.label}</th>)}</tr></thead><tbody>{rows.map((row, index) => <tr key={`${row.rank ?? "row"}-${index}`}>{columns.map((column) => <td key={column.key}>{column.format ? column.format(row[column.key]) : cleanDisplayText(row[column.key]) || "-"}</td>)}</tr>)}</tbody></table></div>;
}

function ResultTabs({ tab, enabled, onChange, className = "" }: { tab: Tab; enabled: boolean; onChange: (tab: Tab) => void; className?: string }) {
  return <nav className={`result-tabs ${className}`} aria-label="Result views">
    {tabs.map((item) => <button key={item} disabled={!enabled} className={tab === item ? "active" : ""} onClick={() => onChange(item)}>{item === "Network" ? <Network size={17} /> : null}{item}</button>)}
  </nav>;
}

function ColumnResizer({
  label,
  value,
  min,
  max,
  direction = 1,
  onResize,
}: {
  label: string;
  value: number;
  min: number;
  max: number;
  direction?: 1 | -1;
  onResize: (next: number) => void;
}) {
  const clamp = (next: number) => Math.min(max, Math.max(min, next));
  const startResize = (event: ReactPointerEvent<HTMLDivElement>) => {
    if (event.button !== 0) return;
    event.preventDefault();
    const startX = event.clientX;
    const startValue = value;
    document.body.classList.add("resizing-columns");
    const move = (moveEvent: globalThis.PointerEvent) => onResize(clamp(startValue + (moveEvent.clientX - startX) * direction));
    const finish = () => {
      document.body.classList.remove("resizing-columns");
      window.removeEventListener("pointermove", move);
      window.removeEventListener("pointerup", finish);
      window.removeEventListener("pointercancel", finish);
    };
    window.addEventListener("pointermove", move);
    window.addEventListener("pointerup", finish);
    window.addEventListener("pointercancel", finish);
  };
  const resizeWithKeyboard = (event: ReactKeyboardEvent<HTMLDivElement>) => {
    if (event.key !== "ArrowLeft" && event.key !== "ArrowRight") return;
    event.preventDefault();
    const delta = (event.key === "ArrowRight" ? 16 : -16) * direction;
    onResize(clamp(value + delta));
  };
  return <div
    className="column-resizer"
    role="separator"
    aria-label={label}
    aria-orientation="vertical"
    aria-valuemin={min}
    aria-valuemax={max}
    aria-valuenow={Math.round(value)}
    tabIndex={0}
    title="Drag or use arrow keys to resize"
    onPointerDown={startResize}
    onKeyDown={resizeWithKeyboard}
  ><span aria-hidden="true" /></div>;
}

export function CytoCartoWorkbench() {
  const [rawText, setRawText] = useState("");
  const [analytes, setAnalytes] = useState<ParsedAnalyte[]>([]);
  const [result, setResult] = useState<CytoCartoResult | null>(null);
  const [tab, setTab] = useState<Tab>("Overview");
  const [inputMode, setInputMode] = useState<InputMode>("paste");
  const [editingIndex, setEditingIndex] = useState<number | null>(null);
  const [numericDrafts, setNumericDrafts] = useState<Record<string, string>>({});
  const [networkCytokines, setNetworkCytokines] = useState<NetworkCytokine[]>([]);
  const [catalogError, setCatalogError] = useState<string | null>(null);
  const [selectedProgram, setSelectedProgram] = useState<string | null>(null);
  const [isRunning, setIsRunning] = useState(false);
  const [isDownloading, setIsDownloading] = useState(false);
  const [helpOpen, setHelpOpen] = useState(false);
  const [analysisPhase, setAnalysisPhase] = useState<"queued" | "running" | null>(null);
  const [error, setError] = useState<string | null>(null);
  const [age, setAge] = useState("");
  const [sex, setSex] = useState("");
  const [race, setRace] = useState("");
  const [ethnicity, setEthnicity] = useState("");
  const [diseaseQuery, setDiseaseQuery] = useState("");
  const [edgeThreshold, setEdgeThreshold] = useState(0.05);
  const [showPositive, setShowPositive] = useState(true);
  const [showNegative, setShowNegative] = useState(true);
  const [showContext, setShowContext] = useState(true);
  const [inputColumnWidth, setInputColumnWidth] = useState(340);
  const [parsedColumnWidth, setParsedColumnWidth] = useState(285);
  const [networkSidebarWidth, setNetworkSidebarWidth] = useState(320);
  const [networkInspectorWidth, setNetworkInspectorWidth] = useState(380);
  const paperUrl = process.env.NEXT_PUBLIC_PAPER_URL;
  const githubUrl = process.env.NEXT_PUBLIC_GITHUB_URL || "https://github.com/poconnel3/CytoCarto";

  useEffect(() => {
    getNetworkCytokines()
      .then((rows) => {
        setNetworkCytokines(rows);
        setCatalogError(null);
      })
      .catch(() => setCatalogError("The DSPIN cytokine list is unavailable. Check the local API."));
  }, []);

  useEffect(() => {
    if (!helpOpen) return;
    const closeOnEscape = (event: globalThis.KeyboardEvent) => {
      if (event.key === "Escape") setHelpOpen(false);
    };
    window.addEventListener("keydown", closeOnEscape);
    return () => window.removeEventListener("keydown", closeOnEscape);
  }, [helpOpen]);

  const clearResults = useCallback(() => {
    setResult(null);
    setTab("Overview");
    setSelectedProgram(null);
  }, []);

  const parseText = useCallback((text: string) => {
    setRawText(text);
    setAnalytes(withEntryIds(parsePanelText(text)));
    setNumericDrafts({});
    setEditingIndex(null);
    setError(null);
    clearResults();
  }, [clearResults]);

  const updateAnalyteName = (index: number, value: string) => {
    setAnalytes((current) => current.map((row, rowIndex) => {
      if (rowIndex !== index) return row;
      return { ...row, raw_name: value };
    }));
    clearResults();
  };

  const updateNumericAnalyte = (index: number, key: NumericAnalyteKey, value: string) => {
    const row = analytes[index];
    if (!row) return;
    setNumericDrafts((current) => ({ ...current, [draftKey(row, index, key)]: value }));
    const parsed = parseDecimalDraft(value);
    setAnalytes((current) => current.map((item, rowIndex) => {
      if (rowIndex !== index) return item;
      if (parsed === null && value !== "") return item;
      const next = {
        ...item,
        [key]: value === "" && key !== "value" ? undefined : parsed ?? item[key],
        qualifier: key === "value" ? "=" as const : item.qualifier,
        flag: undefined,
      };
      return { ...next, flag: derivedAnalyteFlag(next) };
    }));
    clearResults();
  };

  const selectManualAnalyte = (index: number, value: string) => {
    const option = networkCytokines.find((item) => item.display_name === value);
    setAnalytes((current) => current.map((row, rowIndex) => {
      if (rowIndex !== index) return row;
      const next = {
        ...row,
        raw_name: value,
        qualifier: "=" as const,
        reference_low: undefined,
        reference_high: option?.reference_upper_pg_ml ?? undefined,
        flag: undefined,
      };
      return { ...next, flag: derivedAnalyteFlag(next) };
    }));
    const row = analytes[index];
    if (row) {
      setNumericDrafts((current) => ({
        ...current,
        [draftKey(row, index, "reference_low")]: "",
        [draftKey(row, index, "reference_high")]: option?.reference_upper_pg_ml?.toString() ?? "",
      }));
    }
    clearResults();
  };

  const addAnalyte = () => {
    const entryId = crypto.randomUUID();
    setAnalytes((current) => [...current, { raw_name: "", value: 0, qualifier: "=", units: "pg/mL", entry_id: entryId, manual: true }]);
    setNumericDrafts((current) => ({ ...current, [`${entryId}:value`]: "", [`${entryId}:reference_low`]: "", [`${entryId}:reference_high`]: "" }));
    setEditingIndex(analytes.length);
    setInputMode("manual");
    clearResults();
  };

  const loadExample = () => {
    setAge("15");
    setSex("Male");
    setRace("Caucasian");
    setEthnicity("Not Hispanic/Latino");
    setInputMode("paste");
    parseText(EXAMPLE_PANEL_TEXT);
  };

  const uploadFile = async (event: ChangeEvent<HTMLInputElement>) => {
    const file = event.target.files?.[0];
    if (!file) return;
    try {
      const parsed = await parsePanelFile(file);
      setRawText(`Imported ${file.name}`);
      setAnalytes(withEntryIds(parsed));
      setNumericDrafts({});
      setEditingIndex(null);
      setError(null);
      clearResults();
    } catch (uploadError) {
      setError(uploadError instanceof Error ? uploadError.message : "Could not parse the selected file.");
    }
  };

  const analyze = async () => {
    const usableAnalytes = analytes.filter((item, index) => {
      const draft = numericDrafts[draftKey(item, index, "value")];
      return item.raw_name.trim() && Number.isFinite(item.value) && (draft === undefined || parseDecimalDraft(draft) !== null);
    });
    if (!usableAnalytes.length) {
      setError("Add at least one cytokine row before running CytoCarto.");
      return;
    }
    setIsRunning(true);
    setAnalysisPhase("queued");
    setError(null);
    setResult(null);
    try {
      const next = await runAnalysis({
        analytes: usableAnalytes,
        age: age ? Number(age) : undefined,
        sex: sex || undefined,
        race: race || undefined,
        ethnicity: ethnicity || undefined,
      }, setAnalysisPhase);
      setResult(next);
      setTab("Overview");
    } catch (analysisError) {
      setError(analysisError instanceof Error ? analysisError.message : "CytoCarto analysis failed.");
    } finally {
      setIsRunning(false);
      setAnalysisPhase(null);
    }
  };

  const downloadBundle = async () => {
    if (!result || isDownloading) return;
    setIsDownloading(true);
    setError(null);
    try {
      await downloadReport(result, apiPath(result.download_url));
    } catch (downloadError) {
      setError(downloadError instanceof Error ? downloadError.message : "Could not create the CytoCarto report.");
    } finally {
      setIsDownloading(false);
    }
  };

  const programs = asRows(result, "program_scores.tsv");
  const topPrograms = programs.slice(0, 5);
  const cellTypes = asRows(result, "cell_type_enrichment.tsv");
  const genes = asRows(result, "gene_subnetwork.tsv");
  const diseases = asRows(result, "human_disease_hypotheses.tsv");
  const mimics = asRows(result, "human_mimic_hypotheses.tsv");
  const diseaseSearch = diseaseQuery.trim().toLowerCase();
  const matchesDiseaseSearch = (row: TableRow) => !diseaseSearch || Object.values(row).some((value) => String(value ?? "").toLowerCase().includes(diseaseSearch));
  const filteredDiseases = diseases.filter(matchesDiseaseSearch);
  const filteredMimics = mimics.filter(matchesDiseaseSearch);
  const evidenceRows = useMemo(() => {
    const supportingTerms = asRows(result, "disease_supporting_terms.tsv");
    const sourceRows = supportingTerms.length ? supportingTerms : asRows(result, "disease_evidence_edges.tsv");
    const seen = new Set<string>();
    return sourceRows.map((row) => {
      const scores = EVIDENCE_CHANNELS.map(([key, label]) => ({ label, score: Number(row[key] ?? 0) })).filter((item) => Number.isFinite(item.score));
      const strongest = scores.sort((a, b) => b.score - a.score)[0] ?? { label: "-", score: 0 };
      return {
        ...row,
        evidence_library: row.library ?? row.resource ?? "-",
        evidence_term: row.term_label ?? row.disease_name ?? "-",
        evidence_role: cleanDisplayText(row.role ?? row.evidence_kind ?? "-").replaceAll("_", " "),
        evidence_score: strongest.score,
        evidence_support: strongest.label,
        evidence_genes: row.leading_genes ?? row.matched_genes ?? row.patient_gene ?? "-",
      };
    }).filter((row) => {
      const key = `${row.evidence_library}::${row.evidence_term}`;
      if (seen.has(key) || Number(row.evidence_score) <= 0) return false;
      seen.add(key);
      return true;
    }).sort((a, b) => Number(b.evidence_score) - Number(a.evidence_score));
  }, [result]);
  const overviewHypotheses: Array<TableRow & { display_type: "Disease" | "Mimic" }> = [
    ...diseases.slice(0, 3).map((row) => ({ ...row, display_type: "Disease" as const })),
    ...mimics.slice(0, 2).map((row) => ({ ...row, display_type: "Mimic" as const })),
  ];
  const maxCellScore = Math.max(...cellTypes.map((row) => Number(row.score ?? 0)), 0.01);
  const topCytokines = useMemo(() => [...(result?.report.analytes ?? [])]
    .filter((item) => item.elevation_score > 0)
    .sort((a, b) => b.elevation_score - a.elevation_score)
    .slice(0, 6), [result]);
  const networkAnalytes = useMemo(() => {
    const visibleAnalytes = new Set(result?.graph.nodes.filter((node) => node.kind === "analyte").map((node) => node.id.replace(/^analyte:/, "")) ?? []);
    return [...(result?.report.analytes ?? [])]
      .filter((row) => row.elevation_score > 0 && visibleAnalytes.has(row.analyte_id))
      .sort((a, b) => b.elevation_score - a.elevation_score)
      .slice(0, 6);
  }, [result]);
  const selectedProgramValue = selectedProgram ?? (topPrograms[0] ? String(topPrograms[0].program) : null);

  const changeTab = (nextTab: Tab) => {
    setTab(nextTab);
    if (nextTab === "Network" && !selectedProgram && topPrograms[0]) setSelectedProgram(String(topPrograms[0].program));
  };

  const header = <><header className="app-header">
    <div className="brand-lockup" aria-label="CytoCarto, Cytokine Cartography">
      <div className="wordmark">
        <strong>Cyto<span className="wordmark-connector" aria-hidden="true" />Carto</strong>
        <span>Cytokine Cartography</span>
      </div>
      <span className="brand-map" aria-hidden="true"><MapIcon size={56} strokeWidth={1.9} /><ImmuneCellIcon cellType="Dendritic cell" size={24} className="brand-map-cell" /></span>
    </div>
    <div className="header-tools">
      {paperUrl ? <a href={paperUrl} target="_blank" rel="noreferrer"><FileText size={18} /> Paper</a> : null}
      <a className="header-github-link" href={githubUrl} target="_blank" rel="noreferrer" aria-label="View CytoCarto on GitHub" title="View CytoCarto on GitHub"><MarkGithubIcon size={20} aria-hidden="true" /></a>
      <button className="icon-button help-button" onClick={() => setHelpOpen(true)} aria-label="How to interpret CytoCarto output" title="How to interpret CytoCarto output"><CircleHelp size={21} /></button>
      {result ? <button className="button button-secondary report-button" onClick={downloadBundle} disabled={isDownloading}>{isDownloading ? <LoaderCircle className="spin" size={18} /> : <Download size={18} />} {isDownloading ? "Preparing..." : "Download report"}</button> : null}
    </div>
  </header>
  {helpOpen ? <div className="help-backdrop" role="presentation" onMouseDown={(event) => { if (event.target === event.currentTarget) setHelpOpen(false); }}>
    <section className="help-dialog" role="dialog" aria-modal="true" aria-labelledby="interpretation-help-title">
      <div className="help-dialog-header"><div><span>CytoCarto evidence guide</span><h2 id="interpretation-help-title">How to interpret the output</h2></div><button className="icon-button" onClick={() => setHelpOpen(false)} aria-label="Close interpretation help" autoFocus><X size={20} /></button></div>
      <p className="help-intro">Disease and mimic rankings are experimental evidence summaries, not diagnostic probabilities. Interpret them with phenotype, genetics, microbiology, and standard clinical testing.</p>
      <dl className="help-definitions">
        <div><dt>combined_p_value / combined_fdr</dt><dd>Strength and multiple-testing-adjusted significance of direct enrichment across human disease libraries. Smaller values indicate stronger evidence.</dd></div>
        <div><dt>evidence_score</dt><dd><span className="help-formula">−log<sub>10</sub>(combined P value)</span>. Larger values indicate stronger direct disease evidence.</dd></div>
        <div><dt>direct_effect_z</dt><dd>Observed direct-term enrichment relative to degree- and term-size-matched null gene sets.</dd></div>
        <div><dt>dspin_score / chea_score</dt><dd>Agreement after signed propagation through the DSPIN and ChEA-KG directed regulatory networks.</dd></div>
        <div><dt>manifold_score</dt><dd>Similarity within a shared multiview network neighborhood. It is unsigned and does not imply direction or causality.</dd></div>
        <div><dt>matched_genes</dt><dd>Highest-scoring patient genes contributing to the disease or mimic match.</dd></div>
        <div><dt>supporting_sources</dt><dd>Independent direct human disease resources supporting the result.</dd></div>
        <div><dt>counterevidence</dt><dd>Used for Perturb-Seqr evidence when a gene perturbation produces a state opposite to the inferred patient state.</dd></div>
      </dl>
    </section>
  </div> : null}</>;

  if (result && tab === "Network") {
    return <main className="app-shell">
      {header}
      <section className="network-workspace" style={{ "--network-sidebar-width": `${networkSidebarWidth}px`, "--network-inspector-width": `${networkInspectorWidth}px` } as CSSProperties}>
        <aside className="network-sidebar">
          <section>
            <h2>Patient &amp; Panel</h2>
            <dl className="patient-summary">
              <div><dt>Age</dt><dd>{result.report.patient.age ?? "Not provided"}</dd></div>
              <div><dt>Sex</dt><dd>{result.report.patient.sex ?? "Not provided"}</dd></div>
              <div><dt>Race</dt><dd>{result.report.patient.race ?? "Not provided"}</dd></div>
              <div><dt>Ethnicity</dt><dd>{result.report.patient.ethnicity ?? "Not provided"}</dd></div>
              <div><dt>Status</dt><dd>Completed</dd></div>
            </dl>
          </section>
          <section className="network-analytes">
            <h2>Parsed Analytes <span>({networkAnalytes.length})</span></h2>
            <div className="network-analyte-head"><span>Analyte</span><span>Patient value</span></div>
            {networkAnalytes.map((item) => <div className="network-analyte-row" key={item.analyte_id}>
              <span className="elevation-arrow">↑</span><strong>{item.display_name}</strong><span>{item.value_pg_ml.toLocaleString()} pg/mL</span>
            </div>)}
            <p><span className="elevation-arrow">↑</span> Elevated vs reference</p>
          </section>
        </aside>
        <ColumnResizer
          label="Resize patient and panel column"
          value={networkSidebarWidth}
          min={260}
          max={520}
          onResize={(next) => setNetworkSidebarWidth(Math.min(next, Math.max(260, window.innerWidth - networkInspectorWidth - 500)))}
        />
        <section className="network-stage">
          <ResultTabs tab={tab} enabled onChange={changeTab} className="network-tabs" />
          <div className="network-toolbar">
            <label>Edge threshold<select value={edgeThreshold} onChange={(event) => setEdgeThreshold(Number(event.target.value))}><option value="0">All</option><option value="0.05">0.05</option><option value="0.1">0.10</option><option value="0.2">0.20</option></select></label>
            <label><input type="checkbox" checked={showPositive} onChange={(event) => setShowPositive(event.target.checked)} /> Positive <i className="toolbar-edge positive-line" /></label>
            <label><input type="checkbox" checked={showNegative} onChange={(event) => setShowNegative(event.target.checked)} /> Inhibitory <i className="toolbar-edge negative-line" /></label>
            <label><input type="checkbox" checked={showContext} onChange={(event) => setShowContext(event.target.checked)} /> Context <i className="toolbar-edge context-line" /></label>
          </div>
          <NetworkGraph
            graph={{
              nodes: result.graph.nodes,
              edges: result.graph.edges.filter((edge) => {
                if (edge.kind === "cell_context") return showContext;
                if (edge.weight < 0 && !showNegative) return false;
                if (edge.weight >= 0 && !showPositive) return false;
                return edge.kind !== "program_program" || Math.abs(edge.weight) >= edgeThreshold;
              }),
            }}
            selectedProgram={selectedProgramValue}
            onProgramSelect={setSelectedProgram}
          />
          <div className="network-legend"><span><i className="legend-node analyte" />Analyte</span><span><i className="legend-node gene"><Dna size={16} /></i>Gene</span><span><i className="legend-node program" />Program</span><span><ImmuneCellIcon cellType="Monocyte" size={27} />Cell type</span><span><i className="toolbar-edge positive-line" />Positive</span><span><i className="toolbar-edge negative-line" />Inhibitory</span><span><i className="toolbar-edge context-line" />Context</span></div>
        </section>
        <ColumnResizer
          label="Resize gene inspector column"
          value={networkInspectorWidth}
          min={300}
          max={560}
          direction={-1}
          onResize={(next) => setNetworkInspectorWidth(Math.min(next, Math.max(300, window.innerWidth - networkSidebarWidth - 500)))}
        />
        <ProgramGeneDrawer program={selectedProgramValue} result={result} onClose={() => undefined} variant="panel" />
      </section>
      <footer>CytoCarto is an experimental research interpretation tool. It is not a validated diagnostic test.</footer>
    </main>;
  }

  return <main className="app-shell">
    {header}
    <section className="workspace" aria-label="CytoCarto analysis workspace" style={{ "--input-column-width": `${inputColumnWidth}px`, "--parsed-column-width": `${parsedColumnWidth}px` } as CSSProperties}>
      <aside className="input-rail">
        <h1>Patient &amp; Panel</h1>
        <div className="demographics">
          <label>Age<input inputMode="numeric" value={age} onChange={(event) => { setAge(event.target.value); clearResults(); }} placeholder="Optional" /></label>
          <label>Sex<select value={sex} onChange={(event) => { setSex(event.target.value); clearResults(); }}><option value="">Optional</option><option value="Female">Female</option><option value="Male">Male</option></select></label>
          <label>Race<select value={race} onChange={(event) => { setRace(event.target.value); clearResults(); }}><option value="">Optional</option>{RACE_OPTIONS.map(([label, value]) => <option key={value} value={value}>{label}</option>)}</select></label>
          <label>Ethnicity<select value={ethnicity} onChange={(event) => { setEthnicity(event.target.value); clearResults(); }}><option value="">Optional</option>{ETHNICITY_OPTIONS.map((value) => <option key={value} value={value}>{value}</option>)}</select></label>
        </div>
        <div className="rail-divider" />
        <h2>Input mode</h2>
        <div className="input-mode" role="group" aria-label="Cytokine input mode">
          <button className={inputMode === "paste" ? "active" : ""} onClick={() => setInputMode("paste")}><ClipboardPaste size={19} />Paste report</button>
          <button className={inputMode === "upload" ? "active" : ""} onClick={() => setInputMode("upload")}><Upload size={19} />Upload file</button>
          <button className={inputMode === "manual" ? "active" : ""} onClick={() => setInputMode("manual")}><Pencil size={18} />Manual entry</button>
        </div>
        {inputMode === "paste" ? <>
          <label className="field-label" htmlFor="panel-text">Paste cytokine report</label>
          <textarea id="panel-text" value={rawText} onChange={(event) => parseText(event.target.value)} placeholder="Paste a clinical report, CSV, TSV, JSON, or extracted PDF text..." />
        </> : null}
        {inputMode === "upload" ? <label className="upload-dropzone"><Upload size={27} /><strong>Choose a cytokine report</strong><span>TXT, CSV, TSV, JSON, RTF, or text-based PDF</span><input type="file" accept=".txt,.csv,.tsv,.json,.rtf,.pdf" onChange={uploadFile} /></label> : null}
        {inputMode === "manual" ? <div className="manual-entry"><TableProperties size={27} /><p>{catalogError ?? "Choose a DSPIN cytokine, then enter its result and laboratory range in the review table."}</p><button className="button button-secondary" onClick={addAnalyte} disabled={!networkCytokines.length}><Plus size={17} /> Add analyte</button></div> : null}
        <div className="parse-status"><span>{analytes.length ? <><Check size={17} /> Detected {analytes.length} analytes.</> : "No analytes detected."}</span><button onClick={loadExample}>Example</button></div>
        <p className="privacy-note">Do not include names, MRNs, dates of birth, or other identifiers. Raw report text remains in this browser.</p>
        <button className="button button-primary run-button" onClick={analyze} disabled={isRunning || !analytes.length}>{isRunning ? <><LoaderCircle className="spin" size={19} /> {analysisPhase === "queued" ? "Queued..." : "Running analysis..."}</> : <><FlaskConical size={19} /> Analyze</>}</button>
        {error ? <p className="error-message" role="alert">{error}</p> : null}
      </aside>

      <ColumnResizer
        label="Resize patient and panel column"
        value={inputColumnWidth}
        min={280}
        max={520}
        onResize={(next) => setInputColumnWidth(Math.min(next, Math.max(280, window.innerWidth - parsedColumnWidth - 480)))}
      />

      <section className="parsed-panel">
        <div className="parsed-heading"><h2>Parsed Analytes</h2><button className="icon-button" onClick={addAnalyte} aria-label="Add analyte" disabled={!networkCytokines.length}><Plus size={19} /></button></div>
        <div className="parsed-table-wrap">
          <table className="parsed-table">
            <thead><tr><th>Analyte</th><th>Result</th><th>Flag</th><th><span className="sr-only">Edit</span></th></tr></thead>
            <tbody>
              {!analytes.length ? <tr><td colSpan={4} className="parsed-empty">Parsed results will appear here for confirmation.</td></tr> : null}
              {analytes.map((row, index) => {
                const flag = derivedAnalyteFlag(row);
                const valueDraft = numericDrafts[draftKey(row, index, "value")];
                return <Fragment key={row.entry_id ?? `${row.raw_name}-${index}`}>
                <tr>
                  <td>{row.raw_name || "New analyte"}</td>
                  <td>{!row.manual && row.qualifier && row.qualifier !== "=" ? row.qualifier : ""}{valueDraft !== undefined && parseDecimalDraft(valueDraft) === null ? valueDraft || "-" : Number(row.value).toLocaleString()}</td>
                  <td><span className={flag === "H" ? "flag-high" : flag === "L" ? "flag-low" : "flag-normal"}>{flag ?? "—"}</span></td>
                  <td><button className="table-icon" onClick={() => setEditingIndex(editingIndex === index ? null : index)} aria-label={`Edit ${row.raw_name || "analyte"}`}><Pencil size={17} /></button></td>
                </tr>
                {editingIndex === index ? <tr className="edit-row"><td colSpan={4}>
                  <label>Analyte{row.manual ? <select value={row.raw_name} onChange={(event) => selectManualAnalyte(index, event.target.value)}><option value="">Select DSPIN cytokine</option>{networkCytokines.map((cytokine) => <option key={cytokine.stimulus} value={cytokine.display_name}>{cytokine.display_name}</option>)}</select> : <input value={row.raw_name} onChange={(event) => updateAnalyteName(index, event.target.value)} />}</label>
                  <div className="edit-grid"><label>Result<input value={numericDrafts[draftKey(row, index, "value")] ?? row.value} inputMode="decimal" onChange={(event) => updateNumericAnalyte(index, "value", event.target.value)} /></label><label>Lower<input value={numericDrafts[draftKey(row, index, "reference_low")] ?? row.reference_low ?? ""} inputMode="decimal" onChange={(event) => updateNumericAnalyte(index, "reference_low", event.target.value)} /></label><label>Upper<input value={numericDrafts[draftKey(row, index, "reference_high")] ?? row.reference_high ?? ""} inputMode="decimal" onChange={(event) => updateNumericAnalyte(index, "reference_high", event.target.value)} /></label></div>
                  <button className="delete-analyte" onClick={() => { setAnalytes((current) => current.filter((_, rowIndex) => rowIndex !== index)); setEditingIndex(null); clearResults(); }}><Trash2 size={16} /> Remove</button>
                </td></tr> : null}
              </Fragment>;})}
            </tbody>
          </table>
        </div>
        {analytes.length ? <p className="parsed-note">Report reference ranges take precedence over configured defaults.</p> : null}
      </section>

      <ColumnResizer
        label="Resize parsed analytes column"
        value={parsedColumnWidth}
        min={220}
        max={500}
        onResize={(next) => setParsedColumnWidth(Math.min(next, Math.max(220, window.innerWidth - inputColumnWidth - 480)))}
      />

      <section className="results-area">
        <ResultTabs tab={tab} enabled={Boolean(result)} onChange={changeTab} />
        {isRunning ? <div className="analysis-progress"><LoaderCircle className="spin" size={28} /><strong>Running CytoCarto</strong><span>Scoring programs, genes, exact cell types, and human evidence.</span></div> : null}
        {result && tab === "Overview" ? <div className="overview">
          <section className="result-section cytokine-section">
            <div className="section-heading"><h2>Cytokine elevation</h2><span>Positive fold increase vs. upper reference limit</span></div>
            <div className="elevation-chart">
              {topCytokines.map((item) => {
                const ratio = Math.max(item.value_pg_ml / item.reference_upper_pg_ml, 1);
                const position = 2 + Math.min(2, Math.log10(ratio)) * 48;
                return <div className="elevation-row" key={item.analyte_id}>
                  <strong>{item.display_name}</strong>
                  <div className="elevation-track"><i className="reference-line" /><span className="above" style={{ left: "2%", width: `${position - 2}%` }} /><b className="above" style={{ left: `${position}%` }} /></div>
                  <span>{ratio.toFixed(ratio >= 10 ? 1 : 2)}x</span>
                </div>;
              })}
              {topCytokines.length ? <div className="elevation-axis"><span>1x</span><span>10x</span><span>100x+</span></div> : <p className="elevation-empty">No cytokines exceed the upper reference limit.</p>}
            </div>
          </section>
          <section className="result-section programs-section">
            <div className="section-heading"><h2>Top programs <span>(ranked)</span></h2></div>
            <div className="program-strip">
              {topPrograms.map((program, index) => <button key={String(program.program)} onClick={() => setSelectedProgram(String(program.program))} className="program-item">
                <span className="program-rank">{index + 1}</span>
                <div><strong>P{String(program.program_id)}</strong><span>{compactProgramLabel(program)}</span></div>
                <div className="program-score-bar"><i style={{ width: `${Math.max(0, Number(program.score)) * 100}%` }} /></div>
                <b>{Number(program.score).toFixed(3)}</b>
              </button>)}
            </div>
          </section>
          <section className="result-section cells-section">
            <div className="section-heading"><h2>Top cell type drivers</h2><span>{result.report.cell_type_enrichment_mode === "exact_per_cell" ? "Exact 10M-cell enrichment" : "Annotation-inferred"}</span></div>
            <div className="cell-list">
              {cellTypes.map((cell) => <div className="cell-row" key={String(cell.cell_type)}>
                <ImmuneCellIcon cellType={String(cell.cell_type)} size={40} />
                <div className="cell-label"><strong>{String(cell.cell_type)}</strong><span>{Number(cell.fold_enrichment ?? 0).toFixed(2)}x enrichment</span></div>
                <div className="cell-bar"><i style={{ width: `${Number(cell.score ?? 0) / maxCellScore * 100}%` }} /></div>
                <b>{Number(cell.score ?? 0).toFixed(3)}</b>
              </div>)}
            </div>
          </section>
          <section className="result-section disease-section">
            <div className="section-heading"><h2>Disease &amp; mimic evidence</h2><span>Experimental hypothesis layer</span></div>
            <div className="data-table-wrap"><table className="data-table disease-table"><thead><tr><th>Rank</th><th>Condition</th><th>Type</th><th>Evidence</th><th>Matching genes</th><th>Source</th></tr></thead><tbody>
              {overviewHypotheses.map((row, index) => <tr key={`${row.display_type}-${String(row.disease_name)}-${index}`}><td>{index + 1}</td><td>{String(row.disease_name)}</td><td><span className={row.display_type === "Disease" ? "type-disease" : "type-mimic"}>{row.display_type}</span></td><td>{String(row.evidence_tier ?? "-").replaceAll("_", " ")}</td><td>{String(row.matched_genes ?? "-")}</td><td>{String(row.supporting_sources ?? "-")}</td></tr>)}
            </tbody></table></div>
          </section>
        </div> : null}
        {result && tab === "Programs" ? <section className="tab-panel"><div className="section-heading"><h2>All 40 programs</h2><span>Signed DSPIN scores</span></div><div className="data-table-wrap"><table className="data-table"><thead><tr><th>Program</th><th>Score</th><th>Annotation</th><th>Cell context</th></tr></thead><tbody>{programs.map((program) => <tr key={String(program.program)}><td><button className="program-link" onClick={() => setSelectedProgram(String(program.program))}>P{String(program.program_id)}</button></td><td>{Number(program.score).toFixed(3)}</td><td>{cleanDisplayText(program.final_annotation ?? "-")}</td><td>{cleanDisplayText(program.immune_cell_type ?? "-")}</td></tr>)}</tbody></table></div></section> : null}
        {result && tab === "Cell types" ? <section className="tab-panel"><div className="section-heading"><h2>Cell-type enrichment</h2><span>Returned cell types only</span></div><ResultTable rows={cellTypes} columns={[{ key: "rank", label: "Rank" }, { key: "cell_type", label: "Cell type" }, { key: "score", label: "Score", format: (value) => Number(value).toFixed(3) }, { key: "fold_enrichment", label: "Enrichment", format: (value) => `${Number(value).toFixed(2)}x` }, { key: "top_cells", label: "Top-matching cells" }, { key: "mode", label: "Method" }]} /></section> : null}
        {result && tab === "Genes" ? <section className="tab-panel"><div className="section-heading"><h2>Top patient genes</h2><span>Signed network score</span></div><ResultTable rows={genes} columns={[{ key: "gene", label: "Gene" }, { key: "score", label: "Score", format: (value) => Number(value).toFixed(3) }, { key: "direct_cytokine_seed", label: "Direct seed", format: (value) => Number(value).toFixed(3) }, { key: "program_regulator_score", label: "Program regulator", format: (value) => Number(value).toFixed(3) }]} /></section> : null}
        {result && tab === "Disease & mimics" ? <section className="tab-panel split-tables"><label className="table-search"><Search size={19} /><input value={diseaseQuery} onChange={(event) => setDiseaseQuery(event.target.value)} placeholder="Search diseases, mimics, classes, or evidence" aria-label="Search diseases and mimics" /><span>{filteredDiseases.length + filteredMimics.length} results</span></label><div><div className="section-heading"><h2>Named disease evidence</h2><span>Experimental</span></div><ResultTable rows={filteredDiseases} columns={[{ key: "rank", label: "Rank" }, { key: "disease_name", label: "Disease" }, { key: "hypothesis_class", label: "Class" }, { key: "combined_fdr", label: "FDR", format: (value) => Number(value).toExponential(2) }, { key: "evidence_tier", label: "Evidence tier" }]} /></div><div><div className="section-heading"><h2>Acquired mimics</h2><span>Separately ranked</span></div><ResultTable rows={filteredMimics} columns={[{ key: "rank", label: "Rank" }, { key: "disease_name", label: "Mimic" }, { key: "hypothesis_class", label: "Class" }, { key: "combined_fdr", label: "FDR", format: (value) => Number(value).toExponential(2) }, { key: "evidence_tier", label: "Evidence tier" }]} /></div></section> : null}
        {result && tab === "Evidence" ? <section className="tab-panel"><div className="section-heading"><h2>Disease network support</h2><span>Auditable evidence terms</span></div><ResultTable rows={evidenceRows} columns={[{ key: "evidence_library", label: "Library" }, { key: "evidence_term", label: "Term" }, { key: "evidence_role", label: "Role" }, { key: "evidence_score", label: "Strongest score", format: (value) => Number(value).toFixed(3) }, { key: "evidence_support", label: "Strongest channel" }, { key: "evidence_genes", label: "Leading genes" }]} /></section> : null}
      </section>
    </section>
    <footer>CytoCarto is an experimental research interpretation tool. It is not a validated diagnostic test.</footer>
    <ProgramGeneDrawer program={selectedProgram} result={result} onClose={() => setSelectedProgram(null)} />
  </main>;
}
