# Files shipped in inst/extdata

| file | what it is |
|---|---|
| `locked_predictions.csv` | The triage predictions for 50 published signatures, fixed on 2026-09-14 before any external validation (SHA-256 `aea506662efef7149cf336004bb1c0c76349206314fc537f950bdac73aa0d9ca`). |
| `locked_predictions.json` | The model behind them, its two cut-offs and the pre-specified endpoints. |
| `attribution_gastric_GSE183904.csv` | Ten-type cell attribution of each signature in the gastric reference atlas (GSE183904), the features the shipped model was trained on. |
| `attribution_lung_GSE131907.csv` | The same attribution in a lung adenocarcinoma atlas (GSE131907, tumor-tissue cells), computed for the registered test of atlas dependence. |

Read them with `system.file("extdata", "<file>", package = "immunoTriage")`.
