export type ParsedAnalyte = {
  raw_name: string;
  value: number;
  qualifier?: "<" | ">" | "=";
  units: "pg/mL" | "ng/L" | "ng/mL";
  flag?: "H" | "L";
  reference_low?: number;
  reference_high?: number;
  source?: string;
  entry_id?: string;
  manual?: boolean;
};

export type NetworkCytokine = {
  display_name: string;
  stimulus: string;
  reference_upper_pg_ml?: number | null;
};

export type TableRow = Record<string, string | number | boolean | null | undefined>;

export type GraphNode = {
  id: string;
  label: string;
  kind: "analyte" | "gene" | "program" | "cell";
  score?: number;
};

export type GraphEdge = {
  id: string;
  source: string;
  target: string;
  kind: "analyte_gene" | "gene_program" | "program_program" | "cell_context";
  weight: number;
};

export type CytoCartoResult = {
  report: {
    patient: { age?: number | null; sex?: string | null; race?: string | null; ethnicity?: string | null };
    analytes: Array<{
      analyte_id: string;
      display_name: string;
      value_pg_ml: number;
      reference_upper_pg_ml: number;
      elevation_score: number;
      log2_ratio: number;
      stimulus?: string | null;
      stimulus_status: string;
      gene_symbols: string[];
      note?: string;
    }>;
    warnings: string[];
    notes: string[];
    cell_type_enrichment_mode: string;
  };
  tables: Record<string, TableRow[]>;
  download_url: string;
  graph: { nodes: GraphNode[]; edges: GraphEdge[] };
  selected_panel?: string | null;
};
