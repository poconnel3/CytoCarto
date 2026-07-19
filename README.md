# CytoCarto

CytoCarto (Cytokine Cartography) projects a clinical cytokine panel into a fixed human PBMC D-SPIN reference atlas. It provides cytokine normalization, covariate-aware reference responses, signed program and gene-network projection, exact per-cell cell-type enrichment, human disease/mimic evidence retrieval, and human Perturb-Seqr mechanism evidence.

CytoCarto is an experimental research interpretation tool. It is not a validated diagnostic test and does not provide diagnostic probabilities.

## Repository contents

- `web/`: Next.js user interface.
- `web_api/`: FastAPI adapter around the unchanged interpreter.
- `scripts/cytokine_dspin_interpreter.py`: cytokine-to-program, gene, and cell-type projection workflow.
- `scripts/human_kg_disease.py`: human-only disease and mimic evidence layer.
- `scripts/public_evidence.py`: human Perturb-Seqr evidence retrieval and scoring.
- `config/cytokine_panels/`: clinical cytokine-panel aliases and reference defaults.
- `modal_backend.py`: optional Modal ASGI deployment wrapper.

The large D-SPIN reference atlas, per-cell program matrix, cell metadata, public-evidence caches, clinical input files, and generated patient reports are intentionally not included in this repository.

## Local development

1. Place the required D-SPIN reference files in a local `global/` directory and set `DSPIN_GLOBAL_ROOT` to that directory. The interpreter requires the fixed D-SPIN program/gene networks, program representation, program and cell names, gene-to-program matrix, stimulus-response matrices, and aligned `cell_metadata.csv`.
2. Provide the donor cytokine-dictionary workbook through `DSPIN_DONOR_WORKBOOK` when covariate-aware donor weighting is required.
3. Install web dependencies:

   ```bash
   cd web
   npm install
   cd ..
   ```

4. Start the local API and frontend:

   ```bash
   ./scripts/run_cytocarto_local.sh
   ```

5. Open `http://localhost:3000`.

## Tests

```bash
cd web
npm test
npm run lint
npm run build
```

Python interpreter and evidence-layer tests are in `tests/` and API tests are in `web_api/test_app.py`.

## Data and privacy

Raw pasted reports remain in the browser. The API receives normalized cytokine values and optional age, sex, race, and ethnicity. Do not submit names, medical-record numbers, dates of birth, or other direct identifiers.

When using the optional remote Perturb-Seqr refresh path, CytoCarto requires explicit approval before transmitting a patient-derived signed gene list to a third-party public API. Static KG refreshes download public resources only and do not transmit patient data.
