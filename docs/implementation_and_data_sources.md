# Implementace kritérií, klasifikace a datové zdroje ARIANE

## 1. Rozsah aplikace

ARIANE provádí první automatizovaný průchod pravidly ClinGen ENIGMA BRCA1/2 VCEP v1.2 pro geny BRCA1 a BRCA2.

Použité referenční transkripty:

| Gen | RefSeq | Ensembl pro předpočítaná data |
| --- | --- | --- |
| BRCA1 | `NM_007294.4` | `ENST00000357654.9` |
| BRCA2 | `NM_000059.4` | `ENST00000380152.8` |

Automatický výsledek není úplnou expertní klasifikací. Kritéria PS4, PM3, PP1, BS2 a BS4 vyžadují klinická, rodinná nebo literární data a automaticky se nepřidávají. Aplikace pro ně podporuje oddělenou strukturovanou manuální revizi.

Hlavní implementace klasifikace je v `backend/modules/classifier.py`.

## 2. Zpracování vstupu

### 2.1 Normalizace HGVS

Vstup obsahuje gen, `c.` notaci, volitelnou `p.` notaci a u exonové duplikace informaci, zda je uspořádání tandemové nebo neznámé.

Povolené transkripty jsou kontrolovány proti zvolenému genu. Pokud je součástí vstupu například `NM_007294.4:c.303T>G`, prefix transkriptu se ověří a pro interní zpracování se oddělí.

### 2.2 Kontrola referenční alely

Každá SNV prochází před klasifikací kontrolou referenční báze. Kontrola používá lokální verzovaná data, nikoli ruční seznam variant.

Zdroje kontroly:

1. coding SNV snapshot pro kódující pozice,
2. intronická souřadnicová mapa pro podporované intronické pozice.

Příklad:

| Vstup | Výsledek kontroly |
| --- | --- |
| `BRCA1 c.181T>G` | přijato, reference na `c.181` je `T` |
| `BRCA1 c.181A>C` | odmítnuto, reference na `c.181` není `A` |

Neshoda reference vrací HTTP 422 a klasifikace se nespustí. Pokud nainstalovaná data neumožňují referenci ověřit, aplikace postup ukončí místo tichého pokračování.

Implementace je v `backend/modules/reference_validation.py`.

### 2.3 Kontrola proteinového následku

Uživatel musí zadat `c.` i `p.` notaci. Pokud proteinový následek není znám, zadá `p.(?)`. U podporované coding SNV aplikace porovná zadanou `p.` notaci s předpočítaným SNV snapshotem a při rozporu vstup odmítne.

Příklad:

```text
BRCA1 c.303T>G
BRCA1 p.(Tyr101Ter)
variant_type = nonsense
```

Pokud uživatel zadá `p.` notaci, která se liší od následku uloženého pro referenční transkript, aplikace vrátí HTTP 422. Rozpor se neřeší výběrem jedné z hodnot.

### 2.4 Odvození typu varianty

Typ varianty se odvozuje z kombinace `c.` a `p.` notace. Rozlišují se zejména:

- nonsense,
- frameshift,
- missense,
- synonymous,
- intronic,
- canonical splice site,
- in-frame deletion, insertion a delins,
- exon deletion a duplication,
- initiation codon,
- 5' UTR a 3' UTR.

Typ varianty určuje, které větve pravidel lze použít. Nonsense varianta bez proteinového následku se nesmí tiše považovat za missense.

## 3. Automaticky vyhodnocovaná kritéria

Síla kritéria se převádí na body:

| Směr | Síla | Body |
| --- | --- | ---: |
| patogenní | Very Strong | 8 |
| patogenní | Strong | 4 |
| patogenní | Moderate | 2 |
| patogenní | Supporting | 1 |
| benigní | Strong | -4 |
| benigní | Moderate | -2 |
| benigní | Supporting | -1 |

BA1 je samostatné benigní kritérium a vede přímo ke Class 1.

### 3.1 BA1, BS1 a PM2 Supporting

Zdroj: lokální snapshoty gnomAD v2.1.1 exomes non-cancer a gnomAD v3.1.2 genomes non-cancer včetně pokrytí.

Frekvenční metrika se vybírá v pořadí FAF95, popmax a AF podle dostupnosti.

| Kritérium | Podmínka |
| --- | --- |
| BA1 | frekvence nad 0,001 a průměrná hloubka alespoň 20 |
| BS1 Strong | frekvence nad 0,0001 a nejvýše 0,001, hloubka alespoň 20 |
| BS1 Supporting | frekvence nad 0,00002 a nejvýše 0,0001, hloubka alespoň 20 |
| PM2 Supporting | nepřítomnost v požadovaných non-cancer datasetech a průměrná hloubka alespoň 25 |

PM2 se nepoužívá pro indely a exonové CNV. Chybějící záznam bez dostatečného pokrytí není důkazem nepřítomnosti. Fixture nebo neúplná cache nemůže vytvořit frekvenční kritérium.

### 3.2 PVS1 a PM5 PTC

Zdroj: ENIGMA Specifications Table 4 v1.2.

PVS1 se vyhodnocuje pro:

- nonsense varianty,
- frameshift varianty,
- vybrané splice-site varianty,
- exonové delece,
- exonové duplikace.

Table 4 obsahuje pravidla pro jednotlivé exony, kritické C-terminální hranice, splice varianty a exonové přestavby. Výsledná síla může být Very Strong, Strong, Moderate, Supporting nebo N/A.

U splice variant se PVS1 nepřidává pouze podle vzdálenosti od exonu. Varianta musí mít odpovídající pravidlo v Table 4. Větve závislé na RNA jsou označeny k manuální revizi a automatické PVS1 se nepřidá.

U iniciačního kodonu se automatické PVS1 nepoužívá. Aplikace vytvoří doporučení pro strukturovanou revizi podle iniciačního flowchartu.

PM5 PTC se přidává z PTC pravidla stejné Table 4. Síla je určena záznamem pro příslušný exon a typ předčasného terminačního kodonu.

### 3.3 PS3 a BS3

Zdroj: ENIGMA Specifications Table 9 v1.2.

Vyhledávání používá přesný klíč `gene:c_notation`. Automaticky se použijí pouze řádky s přiřazeným PS3 nebo BS3 a podporovanou silou.

Table 9 obsahuje také řádky, ve kterých PS3 ani BS3 nebylo splněno. Tyto řádky zůstávají součástí lossless snapshotu, ale nevytvářejí kritérium.

Pokud má revidovaný řádek Table 9 vlastní hodnotu SpliceAI, používá se tato zmrazená hodnota pro navazující rozhodnutí BP1, BP4 a BP7. Rozdíl proti aktuální cache nebo službě se zobrazí ve varování.

### 3.4 PP4 a BP5

Zdroj: ENIGMA Supplementary Table 7 v1.2.

Kritéria se odvozují z posterior probability a kombinovaného likelihood ratio.

| Kritérium | Supporting | Moderate | Strong | Very Strong |
| --- | ---: | ---: | ---: | ---: |
| PP4 | LR >= 2,08 | LR >= 4,3 | LR >= 18,7 | LR >= 350 |
| BP5 | LR <= 0,48 | LR <= 0,23 | LR <= 0,05 | LR <= 0,00285 |

PP4 je podporováno ve strukturované manuální revizi. Reviewer zadá variantově specifický combined clinical LR, citaci zdroje a souhrn zahrnutých klinických dat včetně kontroly jejich nezávislosti. ARIANE určí sílu výhradně podle uvedených prahů. Sílu PP4 nelze ručně přepsat a neúplný záznam nelze aplikovat. Automatická Module 1 klasifikace PP4 nepřidává.

### 3.5 PS1 na proteinové úrovni

Zdroj referenčních variant: P/LP missense varianty odvozené z ENIGMA Supplementary Table 7.

PS1 vyžaduje:

- missense variantu,
- stejnou aminokyselinovou změnu jako známá P/LP varianta,
- jinou nukleotidovou změnu,
- potvrzené SpliceAI nejvýše 0,1.

Patogenní referenční varianta dává PS1 Strong. Likely Pathogenic referenční varianta dává PS1 Moderate.

PS1 pro stejný splice efekt se automaticky neboduje. Aplikace pouze označí kandidáta pro strukturovanou manuální revizi.

### 3.6 PP3

PP3 Supporting má dvě automatické větve.

SpliceAI větev:

- SpliceAI alespoň 0,2,
- pouze pro povolené typy, například synonymous, missense, in-frame a intronic,
- nepoužívá se jako obecné PVS1 pro nonsense, frameshift, exonové CNV nebo canonical splice-site varianty.

BayesDel_noAF větev:

- missense nebo in-frame varianta v klinicky významné funkční doméně,
- BRCA1: BayesDel_noAF alespoň 0,28,
- BRCA2: BayesDel_noAF alespoň 0,30.

PP3 se nepřičítá současně s aplikovaným PVS1.

### 3.7 BP4

BP4 Supporting vyžaduje potvrzený nízký SpliceAI. Chybějící skóre se nepovažuje za nulové.

Pro missense a in-frame variantu uvnitř funkční domény platí:

- SpliceAI nejvýše 0,1,
- BRCA1 BayesDel_noAF nejvýše 0,15,
- BRCA2 BayesDel_noAF nejvýše 0,18.

Pro synonymous variantu uvnitř domény a pro podporovanou intronickou variantu se používá potvrzené SpliceAI nejvýše 0,1 podle příslušné větve pravidla.

### 3.8 BP7

BP7 Supporting se používá společně s BP4.

U synonymous variant uvnitř funkční domény vyžaduje aplikované BP4 a SpliceAI nejvýše 0,1. U intronických variant se navíc kontroluje, že pozice neleží v konzervovaném donorovém nebo akceptorovém motivu.

Synonymous varianta mimo funkční doménu je řešena přes BP1, nikoli přidáním BP7.

RNA větev BP7 Strong vyžaduje manuální strukturovanou revizi.

### 3.9 BP1

BP1 Strong se používá pro missense, synonymous a in-frame varianty mimo klinicky významnou funkční doménu, pokud je SpliceAI potvrzeně nejvýše 0,1.

Použité domény:

| Gen | Doména | Aminokyselinový rozsah |
| --- | --- | --- |
| BRCA1 | RING | 2 až 101 |
| BRCA1 | coiled-coil | 1391 až 1424 |
| BRCA1 | BRCT | 1650 až 1857 |
| BRCA2 | PALB2 binding | 10 až 40 |
| BRCA2 | DBD | 2481 až 3186 |

Chybějící SpliceAI znamená, že BP1 nelze použít.

## 4. Kritéria vyžadující manuální revizi

Automatická Module 1 klasifikace nepřidává PS4, PM3, PP1, PP4, BS2 a BS4. Tyto kódy závisejí na datech, která nelze bezpečně odvodit pouze z HGVS varianty.

Strukturovaná manuální část dále podporuje:

- PP4 z variantově specifického combined clinical LR,
- PVS1 RNA,
- BP7 RNA,
- PVS1 pro iniciační kodon,
- PS1 pro stejný splice efekt.

Manuálně doplněná kritéria vytvářejí oddělený amended working result. Původní automatická Module 1 klasifikace zůstává zachována.

## 5. Postup klasifikace

Kritéria se vyhodnocují v tomto pořadí:

1. kontrola vstupu, reference a proteinového následku,
2. souřadnice GRCh37 a GRCh38,
3. gnomAD a BA1,
4. Table 9 pro PS3 a BS3,
5. Table 4 pro PVS1 a PM5 PTC,
6. lokální referenční datasety pro navazující variantově specifická pravidla,
7. proteinové PS1,
8. PP3 a BP4,
9. BP7,
10. BP1,
11. klasifikační kombinace.

BA1 ukončí klasifikaci jako Class 1.

Pokud jsou všechna aplikovaná kritéria pouze v jednom směru, používají se kombinace ENIGMA VCEP v1.2 Table 3. Samotný součet bodů v tomto případě neurčuje třídu.

Příklad: PVS1 Very Strong bez dalšího kritéria zůstává Class 3, protože nesplňuje kombinaci pro Likely Pathogenic.

Pokud jsou současně přítomna patogenní i benigní kritéria, používá se bodový systém Tavtigian 2020:

| Součet | Třída |
| ---: | --- |
| 10 a více | Class 5, Pathogenic |
| 6 až 9 | Class 4, Likely Pathogenic |
| -1 až 5 | Class 3, VUS |
| -6 až -2 | Class 2, Likely Benign |
| méně než -6 | Class 1, Benign |

Výsledek s protichůdnými směry obsahuje upozornění na nutnost expertní revize.

## 6. Oficiální ENIGMA tabulky

Oficiální zdroj je ClinGen CSpec GN092, ENIGMA BRCA1/2 VCEP v1.2.0, vydání 2025-01-09:

<https://cspec.genome.network/cspec/ui/svi/doc/GN092/versions>

### 6.1 Table 4

Runtime soubor: `backend/data/enigma_table4.json`

Generátor: `scripts/build_enigma_table4_snapshot.py`

Obsah:

- 493 zdrojových řádků,
- 20 zdrojových sloupců,
- indexy exonových rozsahů,
- PTC pravidla,
- splice pravidla,
- pravidla exonových delecí,
- pravidla exonových duplikací.

### 6.2 Table 9

Runtime soubor: `backend/data/enigma_table9.json`

Generátor: `scripts/build_enigma_table9_snapshot.py`

Obsah:

- 4 731 řádků,
- 14 zdrojových sloupců,
- PS3 a BS3 přiřazení,
- publikované splice výsledky,
- revidované SpliceAI hodnoty,
- 437 revidovaných řádků bez aplikovaného PS3 nebo BS3.

### 6.3 Supplementary Table 7

Runtime soubor: `backend/data/st7_reference_set.json`

Generátor: `scripts/build_enigma_st7_snapshot.py`

Obsah:

- 773 variant,
- 28 zdrojových sloupců,
- prior a posterior probability,
- IARC třída,
- populační a referenční údaje.

### 6.4 Kontrola úplnosti při startu

Table 4, Table 9 a ST7 jsou povinné runtime datasety. `backend/data_validation.py` kontroluje před spuštěním API:

- existenci a čitelnost JSON,
- verzi schématu,
- očekávaný počet řádků a sloupců,
- povinná pole,
- povolené kódy a síly,
- duplicity,
- konzistenci exonových odkazů.

Neúplná nebo poškozená povinná tabulka zastaví start aplikace.

## 7. Předpočítaný coding SNV prostor

### 7.1 Soubor a rozsah

Index:

`data/precomputed/brca_module1_snv_classification_snapshot.index.json`

Metadata:

`data/precomputed/brca_module1_snv_classification_snapshot.metadata.json`

Snapshot obsahuje 47 547 coding SNV pro BRCA1 a BRCA2. Pro každou referenční pozici byly zahrnuty tři možné alternativní báze.

Rozdělení typů ve snapshotu:

| Typ | Počet |
| --- | ---: |
| initiation codon | 18 |
| missense | 35 241 |
| synonymous | 9 891 |
| nonsense | 2 397 |

Každý záznam obsahuje zejména:

- gen a `c.` notaci,
- `p.` notaci,
- typ varianty,
- GRCh37 a GRCh38 souřadnice,
- předpočítané SpliceAI,
- souhrn gnomAD,
- předpočítaná kritéria a třídu z okamžiku vytvoření snapshotu.

### 7.2 Použití za běhu

Runtime používá snapshot pro:

1. kontrola zadané `p.` notace proti následku pro referenční transkript,
2. kontrolu rozporu mezi zadanou a předpočítanou `p.` notací,
3. kontrolu referenční báze coding SNV,
4. lokální převod coding SNV na GRCh37 a GRCh38.

Runtime nepřebírá předpočítanou finální třídu ani seznam kritérií jako hotový výsledek dotazu. Po kontrole vstupu se kritéria znovu vyhodnotí aktuální implementací a aktuálně načtenými runtime datasety.

Toto oddělení umožňuje použít stabilní transkriptový překlad a souřadnice bez zmrazení celé klasifikace na verzi, ve které byl snapshot vytvořen.

### 7.3 Stav a omezení

Metadata označují snapshot jako `snapshot_not_authoritative`. Před klinickým označením za autoritativní dataset je nutná samostatná validace proteinových následků, referenčních alel, souřadnic a reprodukovatelnosti generování.

Tento snapshot pokrývá coding SNV. Coding malé indely pokrývá samostatný snapshot popsaný níže. Exonové CNV a všechny UTR varianty pokryté nejsou.

### 7.4 Normalizovaný snapshot malých indelů

Index:

`data/precomputed/brca_normalized_indel_snapshot.index.json`

Metadata:

`data/precomputed/brca_normalized_indel_snapshot.metadata.json`

Snapshot byl vytvořen z oficiálního BRCA Exchange release 70 ze dne 8. března 2026. Tento release opravil chybná GRCh37 mapování z předchozího release. Zdrojový soubor a výsledný index mají v metadatech SHA-256.

Obsahuje 16 511 jednoznačných malých indelů. Proteinový následek `p.?` zůstává neznámý a nepoužívá se ke kontrole vstupu ani k odvození frameshift nebo in-frame typu:

| Typ | Počet |
| --- | ---: |
| frameshift | 6 898 |
| deletion, proteinový následek neznámý | 4 252 |
| insertion, proteinový následek neznámý | 1 957 |
| duplication, proteinový následek neznámý | 1 897 |
| in-frame deletion | 810 |
| in-frame delins | 269 |
| in-frame duplication | 228 |
| in-frame insertion | 113 |
| delins, proteinový následek neznámý | 87 |

Každý záznam obsahuje vstupní aliasy a kanonickou `c.` notaci, `p.` následek, typ varianty, REF a ALT alely pro GRCh37 a GRCh38, referenční transkript, BRCA Exchange release, zdrojové databáze, CA ID a VRS ID, pokud byly dostupné.

Builder `scripts/build_brca_indel_snapshot.py` čte zdroj streamovaně. Přijímá pouze BRCA1 a BRCA2 na `NM_007294.4` a `NM_000059.4`, malé indely s uvedeným proteinovým následkem, včetně `p.?`, a oběma genomovými mapováními. Konfliktní záznamy se nevkládají. Alias sdílený více záznamy se z aliasového indexu odstraní a zapíše se do metadat. Přesná kanonická notace zůstává dostupná. Aktuální release obsahuje dva takové aliasy. Runtime při startu kontroluje status, počet záznamů a checksum indexu. Chybějící nebo poškozený snapshot zastaví start aplikace.

Za běhu se alias nejprve převede na kanonickou `c.` notaci. Záznam poskytne očekávanou `p.` notaci a lokální souřadnice. Rozpor v `p.` notaci skončí chybou 422. Snapshot neurčuje výslednou klinickou třídu ani automaticky nepřidává kritéria.

## 8. Předpočítaná SpliceAI data

### 8.1 Coding SNV cache

Data:

`data/spliceai/spliceai_brca_snv_reference_cache.json`

Metadata:

`data/spliceai/spliceai_brca_snv_reference_cache.metadata.json`

Obsahuje výsledky pro všech 47 547 variant coding SNV manifestu. Transkriptová politika je `reference_transcript`. Metadata uvádějí 240 kontrol proti veřejnému Broad SpliceAI API a žádný numerický rozdíl mezi úspěšnými odpověďmi.

### 8.2 Intronická cache

Souřadnicová mapa:

`data/coordinates/brca_intronic_snv_coordinates.json`

Metadata:

`data/coordinates/brca_intronic_snv_coordinates.metadata.json`

Mapa obsahuje 13 800 intronických SNV v okně 50 bp od hranic kódujících exonů. Referenční báze pocházejí z UCSC Genome Browser sequence API pro hg19 a hg38. Mapa je navázána na stejné referenční transkripty jako coding snapshot.

SpliceAI výsledky se ukládají do:

- verzované předpočítané snapshoty v `data/spliceai/`, které runtime pouze čte,
- dynamická API cache v `${ARIANE_RUNTIME_CACHE_DIR}/spliceai_api_cache.json`,
- na Railway automaticky v `${RAILWAY_VOLUME_MOUNT_PATH}/ariane-runtime-cache/spliceai_api_cache.json`.

Zápis dynamické cache je atomický. Bez nakonfigurovaného runtime adresáře nebo
Railway volume se při lokálním vývoji používá původní `data/spliceai/`.

Intronický předpočítaný snapshot:

`data/spliceai/spliceai_brca_intronic_snv_reference_cache.json`

Cache je použitelná pouze s kompletními metadaty a úspěšnou kontrolou počtu záznamů a checksumu. Rozpracovaná cache se nenačte jako platný zdroj.

### 8.3 Priorita zdrojů

Pro podporované varianty se používá lokální předpočítaná cache. Stav zdroje a důvod selhání se zaznamenávají. Chybějící skóre se nepřevádí na nulu a nemůže vytvořit BP1, BP4 ani BP7.

## 9. Souřadnice

Coding SNV používají přednostně GRCh37 a GRCh38 ze SNV snapshotu. Intronické SNV používají rozšířenou lokální souřadnicovou mapu.

Pro varianty mimo tyto mapy existuje resolver VariantValidator a sekundární resolver Mutalyzer. Stav resolveru, neúplná odpověď a použitý zdroj se zapisují do diagnostiky. Nejednoznačný výsledek se nesmí převést na první nalezené ID.

Souřadnice jsou verzovaná vlastnost kombinace referenčního transkriptu a genomového sestavení. Nemění se při každém dotazu, ale musí se znovu validovat při změně transkriptu, sestavení nebo normalizační politiky.

## 10. gnomAD data

Runtime soubory jsou v `backend/data/gnomad/`.

Používané datasety:

- gnomAD v2.1.1 exomes non-cancer na GRCh37,
- gnomAD v3.1.2 genomes non-cancer na GRCh38,
- samostatná per-position coverage cache.

V3 snapshot se sestavuje skriptem `scripts/build_gnomad_v3_brca_snapshot.py` z oficiálního regionálního VCF gnomAD v3.1.2 a z gnomAD browser coverage API. Stávající v2.1 záznamy se při sestavení zachovávají.

PM2 vyžaduje prokázanou nepřítomnost a dostatečné pokrytí. Samotná absence varianty v JSON není dostačující.

## 11. BayesDel a AlphaMissense

BayesDel_noAF a AlphaMissense se získávají jedním dotazem nad genomovou variantou přes MyVariant.info a ukládají se do lokální cache.

BayesDel se používá pouze v rozhodovacích větvích PP3 a BP4 popsaných výše. AlphaMissense se vrací jako doplňující anotace a samo o sobě nevytváří samostatné ENIGMA kritérium.

Selhání služby, chybějící GRCh37 souřadnice nebo nenalezená anotace mají odlišné stavové kódy. Důvod se zobrazí v diagnostice. Chybějící BayesDel nemůže být nahrazen předpokládanou hodnotou.

## 12. ClinVar a ClinGen

ClinVar a ClinGen ERepo se používají pro externí srovnání a auditní kontext. Jejich klasifikace se automaticky nepřičítá jako ACMG nebo ENIGMA kritérium.

Pokud vyhledávání vrátí více kandidátů bez jednoznačné přesné shody, stav je `ambiguous`. Aplikace nevybere první ID.

## 13. Degradované datové zdroje a fail-closed chování

Povinné ENIGMA tabulky zastaví start aplikace, pokud jsou neúplné nebo nečitelné.

U ostatních zdrojů se degradace ukládá do centrálního registru a zobrazuje v klasifikačním výsledku a `/api/health`. Hlášení obsahuje komponentu a konkrétní důvod. Absolutní serverové cesty se pro uživatele zkracují na tvar `…ariane/...`.

Bez potvrzeného vstupu nebo bez dat požadovaných konkrétním pravidlem se dané kritérium nepoužije. Aplikace nesmí nahrazovat chybějící data nulou, fixture hodnotou, prvním nalezeným ID nebo ručním fallback slovníkem.

## 14. Audit a reprodukovatelnost

Každý předpočítaný dataset má být doprovázen metadaty obsahujícími alespoň:

- název a verzi datasetu,
- datum vytvoření,
- zdroj,
- referenční transkripty a genomová sestavení,
- počet záznamů,
- checksum,
- stav úplnosti,
- verzi nebo commit generátoru, pokud je k dispozici.

Změna oficiální ENIGMA verze, transkriptu, genomového sestavení, predikčního modelu nebo klasifikační logiky vyžaduje nové sestavení dotčených snapshotů a opakování regresních testů.

Související dokumenty:

- `docs/enigma_source_data_audit.md`,
- `docs/manual_evidence_review.md`,
- `docs/splice_ps1_reference_set.md`,
- `docs/vus_explanation_and_golden_cases.md`.
