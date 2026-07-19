"""Modal deployment wrapper for the existing CytoCarto FastAPI adapter."""

from __future__ import annotations

import os
from typing import Any

import modal


APP_NAME = "cytocarto-api"
VOLUME_NAME = "cytocarto-dspin-global"
VOLUME_ROOT = "/mnt/cytocarto"

volume = modal.Volume.from_name(VOLUME_NAME)
image = (
    modal.Image.debian_slim(python_version="3.11")
    .pip_install_from_requirements("web_api/requirements.txt")
    .env(
        {
            "DSPIN_GLOBAL_ROOT": f"{VOLUME_ROOT}/global",
            "DSPIN_EVIDENCE_CACHE": f"{VOLUME_ROOT}/cache/public_evidence_cache",
            "DSPIN_KG_CACHE": f"{VOLUME_ROOT}/cache/human_kg_cache",
            "DSPIN_DONOR_WORKBOOK": f"{VOLUME_ROOT}/inputs/cytokine_dictionary.xlsx",
            "CYTOCARTO_PANEL_CONFIG_DIR": f"{VOLUME_ROOT}/config/cytokine_panels",
        }
    )
    .add_local_dir("scripts", remote_path="/root/scripts", copy=True)
    .add_local_dir("web_api", remote_path="/root/web_api", copy=True)
    .add_local_dir("config/cytokine_panels", remote_path="/root/web_api/panel_configs", copy=True)
    .add_local_dir("config", remote_path="/root/config", copy=True)
)

app = modal.App(APP_NAME)

worker_options = {
    "image": image,
    "volumes": {VOLUME_ROOT: volume.with_mount_options(read_only=True)},
    "memory": 16_384,
    "timeout": 900,
    "max_containers": 2,
    "scaledown_window": 300,
}


@app.function(**worker_options)
def analyze_job(payload: dict[str, Any]) -> dict[str, Any]:
    from web_api.app import AnalysisRequest, run_request

    return run_request(AnalysisRequest.model_validate(payload))


@app.function(**worker_options)
@modal.asgi_app()
def fastapi_app():
    os.environ["DSPIN_GLOBAL_ROOT"] = f"{VOLUME_ROOT}/global"
    os.environ["DSPIN_EVIDENCE_CACHE"] = f"{VOLUME_ROOT}/cache/public_evidence_cache"
    os.environ["DSPIN_KG_CACHE"] = f"{VOLUME_ROOT}/cache/human_kg_cache"
    os.environ["DSPIN_DONOR_WORKBOOK"] = f"{VOLUME_ROOT}/inputs/cytokine_dictionary.xlsx"
    os.environ["CYTOCARTO_PANEL_CONFIG_DIR"] = f"{VOLUME_ROOT}/config/cytokine_panels"
    os.environ["CYTOCARTO_ALLOWED_ORIGINS"] = ",".join(
        [
            "https://cytocarto.vercel.app",
            "http://localhost:3000",
            "http://127.0.0.1:3000",
        ]
    )

    from web_api import app as api_module

    def submitter(payload: dict[str, Any]) -> str:
        return analyze_job.spawn(payload).object_id

    def status_provider(job_id: str) -> dict[str, Any]:
        call = modal.FunctionCall.from_id(job_id)
        try:
            result = call.get(timeout=0)
        except TimeoutError:
            return {"job_id": job_id, "status": "running"}
        except Exception as error:
            return {"job_id": job_id, "status": "failed", "detail": f"CytoCarto analysis failed: {error}"}
        return {"job_id": job_id, "status": "complete", "result": result}

    api_module.configure_job_backend(submitter, status_provider)

    return api_module.app
