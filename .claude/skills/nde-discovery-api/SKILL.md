---
name: nde-discovery-api
description: Query the NIAID Data Ecosystem Discovery API to search and retrieve biomedical datasets, computational tools, clinical studies, and related resources from NIAID-affiliated data repositories.
license: MIT
compatibility: opencode
metadata:
  base_url: https://api.data.niaid.nih.gov/v1
  staging_url: https://api-staging.data.niaid.nih.gov/v1
  docs: https://data.niaid.nih.gov
---

## What I do

I provide instructions for querying the NIAID Data Ecosystem Discovery API at `https://api.data.niaid.nih.gov/v1`. This API indexes biomedical resources including datasets, computational tools, and publications from NIAID-affiliated repositories.

**Scale (verified 2026-09-01, prod build `20260830`):** **14.2M** records across **48** source repositories.

| `@type` | prod | staging |
|---|--:|--:|
| `Sample` | 8,729,088 | ✓ |
| `Dataset` | 5,442,683 | ✓ |
| `ComputationalTool` | 33,272 | ✓ |
| `ResourceCatalog` | 43 | ✓ |
| `Inference` | **0 — not in prod** | **10,399,895** |
| `DataCollection` | **0 — not in prod** | 749,788 |

> 🚨 **`Inference` and `DataCollection` are STAGING-ONLY.** They do not exist in production. An
> `@type:Inference` query against `api.data.niaid.nih.gov` returns **0 hits with no error** — it
> fails silently. You must use `https://api-staging.data.niaid.nih.gov/v1` for these.
> (`Inference` = Gene Expression Atlas 10,395,623 + ImmuneSpace 4,272. Added by nde-crawlers
> PR #368, merged 2026-07-10; not yet promoted to prod as of 2026-09-01. Re-check before relying
> on staging — and note staging may be rebuilt without warning.)

> ⚠️ **`Sample` is easy to miss.** Most GEO content is indexed at **GSM (sample) granularity** as
> well as GSE — GEO alone is 8.8M records, 62% of the prod index. Always pin `@type` or samples
> will dominate your results.

## API Base URL

```
https://api.data.niaid.nih.gov/v1
```

> **Note:** HTTP redirects to HTTPS. Always use `https://` to avoid an extra round-trip.

Staging (for testing):
```
https://api-staging.data.niaid.nih.gov/v1
```

---

## Endpoints

### 1. GET /query — Search resources

Full-text and field-specific search across all indexed resources.

**URL:** `GET https://api.data.niaid.nih.gov/v1/query`

**Key parameters:**

| Parameter | Type | Description |
|---|---|---|
| `q` | string | Query string. Free-text or Lucene syntax. Examples: `"SARS-CoV-2"`, `"@type:Dataset"`, `"infectiousAgent.name:HIV"` |
| `fields` | string | Comma-separated fields to return. Use `fields=all` for everything. E.g. `name,description,url,datePublished` |
| `aggs` | string | Fields to aggregate/facet on. E.g. `includedInDataCatalog.name,infectiousAgent.name` |
| `facet_size` | integer | Number of facet buckets (1–1000, default 10) |
| `size` | integer | Number of hits to return (default 10, max 1000 with `fetch_all`) |
| `from` | integer | Offset for pagination (default 0). **Hard cap 10,000** — `from=10001` returns HTTP 400. Beyond that use `fetch_all=true` + `scroll_id` (500 hits/page) |
| `sort` | string | Comma-separated sort fields. Prefix `-` for descending. E.g. `-datePublished` |
| `extra_filter` | string | Additional Lucene filter. E.g. `(healthCondition.name:("Asthma"))+AND+(includedInDataCatalog.name:("Vivli"))` |
| `fetch_all` | boolean | Return up to 1000 hits when `true` |
| `scroll_id` | string | For paginating through large result sets |
| `hist` | string | Date field to histogram by. E.g. `datePublished` |
| `hist_interval` | string | Histogram interval: `day`, `week`, `month`, `quarter`, `year` |
| `explain` | boolean | Include scoring explanation |

**Example requests:**

```bash
# Free-text search
curl "https://api.data.niaid.nih.gov/v1/query?q=COVID-19&fields=name,description,url&size=5"

# Search for datasets only
curl "https://api.data.niaid.nih.gov/v1/query?q=%40type:Dataset%20AND%20SARS-CoV-2&fields=name,url,includedInDataCatalog.name&size=10"

# Faceted search by repository
curl "https://api.data.niaid.nih.gov/v1/query?q=influenza&aggs=includedInDataCatalog.name&facet_size=20"

# Filter by infectious agent and repository
curl "https://api.data.niaid.nih.gov/v1/query?q=*&extra_filter=(infectiousAgent.name:(%22HIV%22))+AND+(includedInDataCatalog.name:(%22ImmPort%22))"

# Paginate results
curl "https://api.data.niaid.nih.gov/v1/query?q=malaria&size=10&from=10"

# Sort by date descending
curl "https://api.data.niaid.nih.gov/v1/query?q=tuberculosis&sort=-datePublished&fields=name,datePublished,url"

# Histogram by year published
curl "https://api.data.niaid.nih.gov/v1/query?q=Ebola&hist=datePublished&hist_interval=year"
```

**Response structure:**
```json
{
  "took": 42,
  "total": 1234,
  "hits": [
    {
      "_id": "ZENODO_123456",
      "@type": "Dataset",
      "name": "...",
      "description": "...",
      "url": "...",
      "datePublished": "2023-01-15",
      "includedInDataCatalog": [{ "name": "Zenodo" }],
      "infectiousAgent": [{ "name": "Severe acute respiratory syndrome coronavirus 2", "alternateName": ["SARS-CoV-2", "2019-nCoV"] }]
    }
  ],
  "facets": {
    "includedInDataCatalog.name": {
      "_type": "terms",
      "terms": [{ "term": "Zenodo", "count": 425114 }],
      "other": 82658,
      "missing": 317,
      "total": 4390834
    }
  }
}
```

> **Important:** `includedInDataCatalog` is returned as an **array**. Facet buckets use `term`/`count` keys (not `key`/`doc_count`).

---

### 2. POST /query — Batch query

**URL:** `POST https://api.data.niaid.nih.gov/v1/query`

**Query parameters:** same `q`, `scopes`, `from`, `sort` as GET.

**Request body (JSON):**
```json
{
  "q": "HIV",
  "scopes": "_id",
  "from": 0,
  "sort": "-datePublished"
}
```

---

### 3. GET /metadata — API metadata

Returns information about the indexed data sources and statistics.

**URL:** `GET https://api.data.niaid.nih.gov/v1/metadata`

```bash
curl "https://api.data.niaid.nih.gov/v1/metadata"
```

---

### 4. GET /metadata/fields — Available fields

Returns all searchable/returnable metadata fields with descriptions.

**URL:** `GET https://api.data.niaid.nih.gov/v1/metadata/fields`

**Parameters:**

| Parameter | Type | Description |
|---|---|---|
| `search` | string | Filter fields by keyword. E.g. `search=author` |
| `prefix` | string | Filter fields by prefix. E.g. `prefix=study` |
| `format` | string | Output format: `json` (default) or `yaml` |
| `raw` | boolean | Return plain text if `true` |

```bash
# List all fields
curl "https://api.data.niaid.nih.gov/v1/metadata/fields"

# Search for author-related fields
curl "https://api.data.niaid.nih.gov/v1/metadata/fields?search=author"

# Fields with 'infectious' prefix (returns 17 fields including infectiousAgent.name, .url, etc.)
curl "https://api.data.niaid.nih.gov/v1/metadata/fields?search=infectious"
```

---

## Key metadata fields for querying

| Field | Description | Example values |
|---|---|---|
| `@type` | Resource type | `Sample` (8.7M), `Dataset` (5.4M), `ComputationalTool` (33K), `ResourceCatalog` (43); staging-only: `Inference` (10.4M), `DataCollection` (750K) — **no ClinicalTrial type** |
| `name` | Resource title | |
| `description` | Full description | |
| `url` | Resource URL | |
| `_id` | Unique ID (source-prefixed) | `ZENODO_188418`, `OMICSDI_S-EPMC6271715` |
| `datePublished` | Publication date | `2023-01-15` |
| `includedInDataCatalog.name` | Source repository | `Figshare` (1.7M), `NCBI BioProject` (1M), `Zenodo` (425K), `NCBI SRA` (391K), `NCBI GEO` (264K), `Mendeley`, `Protein Data Bank`, `Harvard Dataverse`, `bio.tools`, `Vivli`, `HuBMAP` |
| `infectiousAgent.name` | Pathogen (canonical, lowercase) | `severe acute respiratory syndrome coronavirus 2`, `human immunodeficiency virus`, `influenza a virus` — use `aggs=infectiousAgent.name` to discover exact values |
| `infectiousDisease.name` | Disease | `COVID-19`, `Tuberculosis` |
| `healthCondition` | Health condition | `Asthma`, `Cancer` |
| `species.name` | Species studied | `Homo sapiens`, `Mus musculus` |
| `author.name` | Author name | |
| `keywords` | Keywords/tags | |
| `measurementTechnique.name` | Measurement method | `RNA-seq`, `flow cytometry` |
| `funding.funder.name` | Funding organization | `NIH`, `NIAID` |
| `studyDesign.studyType` | Study type | `Observational Study`, `Interventional` |
| `studyStatus.status` | Trial status | `Recruiting`, `Completed` |
| `doi` | DOI | |
| `pmid` | PubMed ID | |
| `hasDownload` | Availability of download | `true`, `false` |


A complete list of available metadata fields can be found in this [API endpoint](https://api.data.niaid.nih.gov/v1/metadata/fields).
---

## Query syntax

The `q` parameter supports Lucene query syntax:

```
# Match any field
q=SARS-CoV-2

# Field-specific search
q=name:COVID

# Boolean operators
q=influenza AND vaccine

# Type filter
q=@type:Dataset

# Range query
q=datePublished:[2020 TO 2023]

# Wildcard
q=influ*

# Phrase match
q="reference genome"

# Combined
q=@type:ComputationalTool AND measurementTechnique.name:RNA-seq
```

---

## Common use cases

### Discover exact infectiousAgent values before filtering
```bash
# Always do this first to find the canonical name used in the index
curl "https://api.data.niaid.nih.gov/v1/query?q=SARS-CoV-2&aggs=infectiousAgent.name&facet_size=10&size=0"
# Returns terms like: "severe acute respiratory syndrome coronavirus 2"
```

### Find datasets about a pathogen (use canonical name)
```bash
curl "https://api.data.niaid.nih.gov/v1/query?q=%40type:Dataset%20AND%20infectiousAgent.name:%22severe%20acute%20respiratory%20syndrome%20coronavirus%202%22&fields=name,url,datePublished,includedInDataCatalog.name&size=10"
```

### Explore available repositories
```bash
curl "https://api.data.niaid.nih.gov/v1/query?q=*&aggs=includedInDataCatalog.name&facet_size=20&size=0"
# Note: facet buckets use "term" and "count" keys (not "key"/"doc_count")
```

### Find computational tools
```bash
curl "https://api.data.niaid.nih.gov/v1/query?q=%40type:ComputationalTool%20AND%20RNA-seq&fields=name,description,url,codeRepository&size=10"
```

### Find HIV datasets sorted by most recent
```bash
curl "https://api.data.niaid.nih.gov/v1/query?q=%40type:Dataset%20AND%20HIV&fields=name,datePublished,url,includedInDataCatalog.name&sort=-datePublished&size=10"
```

### Get resource by ID
```bash
curl "https://api.data.niaid.nih.gov/v1/query?q=_id:ZENODO_188418"
```

### Publication date histogram
```bash
curl "https://api.data.niaid.nih.gov/v1/query?q=tuberculosis&hist=datePublished&hist_interval=year&size=0"
```

---

## `@type:Inference` — precomputed differential-expression contrasts (STAGING ONLY)

One record per **DE contrast** (gene × comparison), not per dataset. 10.4M records, almost all from
Gene Expression Atlas. This is the closest thing NDE has to a machine-readable
`perturbation → gene` assertion.

```bash
B=https://api-staging.data.niaid.nih.gov/v1     # NOT prod — see the warning at the top
curl -sG "$B/query" --data-urlencode 'q=@type:Inference AND observationAbout.name:TLR7' \
                    --data-urlencode 'size=0'
```

**Record shape** (`gxa_e_mtab_7745_g2_g1_ensg00000142319`):

| Field | Meaning | Example |
|---|---|---|
| `observationAbout` | the **gene** (structured) | `{"name":"SLC6A3","identifier":"ENSG00000142319"}` |
| `value` + `unitText` | effect size | `12.3`, `"Log2 fold change"` |
| `marginOfError.value` | adjusted p-value | `0.0` |
| `measurementQualifier` | the contrast, free text | `'dexamethasone; 1 micromolar' vs 'vehicle'` |
| `variableMeasured` | **test** arm (`StatisticalVariable`) | `.value`, `.description` |
| `measurementDenominator` | **reference** arm | same shape |
| `semanticMapping` | Biolink-typed triple | `triplePredicate.identifier` = `…DirectionQualifierEnum#upregulated`; `tripleSubjectQualifier` = `aspect: abundance`, `direction: increase` |
| `healthCondition` | MONDO-mapped (disease contrasts) | `{"identifier":"0005005","inDefinedTermSet":"MONDO"}` |
| `subjectOf` | parent experiment + **GSM accessions** | `{"identifier":["E-MTAB-7745","GSM…"]}` |

### 🚨 Matching a drug/compound: two mandatory filters

There is **no chemical field**. The compound appears only as free text inside `variableMeasured` /
`measurementDenominator` / `measurementQualifier`. Naive matching produces wrong answers:

```bash
# WRONG — counts contrasts where the drug is in BOTH arms
q=@type:Inference AND dexamethasone

# RIGHT — drug in the test arm, absent from the reference arm, human only
q=@type:Inference AND species.identifier:9606
  AND variableMeasured.value:dexamethasone
  AND NOT measurementDenominator.value:dexamethasone
```

Why it matters: `gxa_e_geod_54049_g1_g5_at1g03790` contrasts
`'pBeaconRFP_GR::bZIP1; …10 uM dexamethasone…'` vs `'empty vector; …10 uM dexamethasone…'` — in
*Arabidopsis*, where dexamethasone is a **GR inducer**, not the variable. The real contrast is
genotype. The `NOT measurementDenominator` clause excludes **7,590** such dexamethasone contrasts;
adding `species.identifier:9606` leaves **23,314** genuine human drug-vs-vehicle contrasts.

Expect the same trap for other tool compounds (tamoxifen/Cre, doxycycline/Tet-on).

---

## When to use me

Use this skill when the user wants to:
- Search NIAID biomedical data resources for datasets, tools, trials, or protocols
- Explore available repositories in the NIAID Data Ecosystem
- Filter resources by infectious agent, disease, species, or data type
- Retrieve metadata fields or API schema information
- Build queries against the NDE Discovery API

Always prefer `curl` or Python `requests` for actual API calls. Use `fields=` to limit response size when only specific attributes are needed. Use `aggs=` to explore the distribution of values in any field before constructing precise queries.
