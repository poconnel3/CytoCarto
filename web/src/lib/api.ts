import type { CytoCartoResult, NetworkCytokine, ParsedAnalyte } from "@/types/cytocarto";

const apiUrl = process.env.NEXT_PUBLIC_CYTOCARTO_API_URL ?? "http://127.0.0.1:8000";

export function apiPath(path: string): string {
  if (/^https?:\/\//i.test(path)) return path;
  return `${apiUrl.replace(/\/$/, "")}/${path.replace(/^\//, "")}`;
}

async function fetchJson<T>(path: string, init?: RequestInit): Promise<T> {
  const response = await fetch(`${apiUrl}${path}`, init);
  if (!response.ok) {
    const payload = (await response.json().catch(() => null)) as { detail?: string } | null;
    throw new Error(payload?.detail ?? `Request failed (${response.status}).`);
  }
  return response.json() as Promise<T>;
}

export function getDemo(): Promise<CytoCartoResult> {
  return fetchJson<CytoCartoResult>("/v1/demo");
}

type AnalysisPayload = { analytes: ParsedAnalyte[]; age?: number; sex?: string; race?: string; ethnicity?: string };
type AnalysisJobResponse = { job_id: string; status: "queued" };
type AnalysisJobStatus = { job_id: string; status: "queued" | "running" | "complete" | "failed"; result?: CytoCartoResult; detail?: string };
type AnalysisStatusCallback = (status: "queued" | "running") => void;

export function normalizedAnalytesForApi(analytes: ParsedAnalyte[]): Array<Omit<ParsedAnalyte, "source">> {
  return analytes.map((analyte) => ({
    raw_name: analyte.raw_name,
    value: analyte.value,
    qualifier: analyte.qualifier,
    units: analyte.units,
    flag: analyte.flag,
    reference_low: analyte.reference_low,
    reference_high: analyte.reference_high,
  }));
}

export async function runAnalysis(payload: AnalysisPayload, onStatus?: AnalysisStatusCallback): Promise<CytoCartoResult> {
  const requestBody = JSON.stringify({ ...payload, analytes: normalizedAnalytesForApi(payload.analytes) });
  const job = await fetchJson<AnalysisJobResponse>("/v1/jobs", {
    method: "POST",
    headers: { "Content-Type": "application/json" },
    body: requestBody,
  });
  onStatus?.(job.status);

  const deadline = Date.now() + 20 * 60 * 1000;
  while (Date.now() < deadline) {
    const response = await fetch(`${apiUrl}/v1/jobs/${encodeURIComponent(job.job_id)}`);
    const payload = (await response.json().catch(() => null)) as AnalysisJobStatus | { detail?: string } | null;
    if (response.status === 202) {
      onStatus?.((payload as AnalysisJobStatus | null)?.status === "queued" ? "queued" : "running");
      await new Promise((resolve) => window.setTimeout(resolve, 2500));
      continue;
    }
    if (!response.ok) {
      throw new Error((payload as { detail?: string } | null)?.detail ?? `Request failed (${response.status}).`);
    }
    const status = payload as AnalysisJobStatus;
    if (status.status === "complete" && status.result) return status.result;
    if (status.status === "failed") throw new Error(status.detail ?? "CytoCarto analysis failed.");
    onStatus?.(status.status === "queued" ? "queued" : "running");
    await new Promise((resolve) => window.setTimeout(resolve, 2500));
  }
  throw new Error("The analysis is still running. Please try again shortly.");
}

export function getProgramGenes(program: string): Promise<Array<{ gene: string; program_weight: number; direction: string }>> {
  return fetchJson<{ genes: Array<{ gene: string; program_weight: number; direction: string }> }>(`/v1/program-genes?program=${encodeURIComponent(program)}`).then((payload) => payload.genes);
}

export function getNetworkCytokines(): Promise<NetworkCytokine[]> {
  return fetchJson<{ cytokines: NetworkCytokine[] }>("/v1/network-cytokines").then((payload) => payload.cytokines);
}
