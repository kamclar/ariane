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

### 2.1 Normalizace vstupu a HGVS

Uživatel zadává gen a jednu variantu. Vstupní normalizační vrstva
`backend/modules/variant_input.py` přijímá:

- referenční transkriptovou notaci, například `c.303T>G`,
- notaci s accession prefixem, například `NM_007294.4:c.303T>G`,
- kombinovanou starší formu s `p.` následkem,
- genomickou variantu ve tvaru `chr17:43099813:C>T`, `17:43099813 C>T`
  nebo `17-43099813-C-T`.

Pro genomickou variantu je povinná sestava `GRCh37` nebo `GRCh38`. Aplikace ji
neodhaduje. Genomický vstup se hledá v obousměrném indexu vytvořeném z
verzovaného coding SNV snapshotu a normalizovaného indel snapshotu. Jedna
jednoznačná shoda vrátí kanonickou `c.` notaci, třípísmennou `p.` notaci,
referenční transkript a zdroj normalizace. Nulová nebo víceznačná shoda ukončí
požadavek před klasifikací.

Povolené transkripty jsou kontrolovány proti zvolenému genu. Aktuální registr
obsahuje `BRCA1: NM_007294.4` a `BRCA2: NM_000059.4`. Přidání dalšího genu
vyžaduje deklarovaný referenční transkript a validovaný normalizační snapshot,
nikoli ruční slovník jednotlivých variant.

Samostatné API `POST /api/normalize` používá stejnou vrstvu jako klasifikační
endpoint. Vrací zadanou notaci, kanonickou `c.` a `p.` notaci, transkript,
sestavu a zdroj. Klasifikace proto normalizaci nemůže obejít.

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

### 2.3 Odvození a kontrola proteinového následku

U coding SNV a normalizovaných indelů aplikace odvodí `p.` notaci z lokálního
snapshotu pro referenční transkript. U intronického nebo UTR vstupu bez
deterministického proteinového následku použije `p.(?)`. Uživatel tedy nemusí
`p.` notaci opisovat.

Veřejný výsledek u `p.(?)` zobrazuje samostatné vysvětlení. Rozlišuje, že
varianta leží mimo překládanou kódující sekvenci nebo může ovlivnit splicing a
že proteinový následek nelze určit pouze z DNA notace. Upozornění uvádí potřebu
transkriptové, RNA, breakpointové nebo strukturální evidence podle typu
varianty. Otazník se proto nezobrazuje bez kontextu jako domnělá chyba parseru.

Starší API a batch vstup mohou `p.` notaci stále dodat. V takovém případě se
porovná s lokálním referenčním následkem. Rozpor vrátí HTTP 422 a klasifikace se
nespustí. Pokud aplikace pro coding variantu nemá validovaný následek a uživatel
jej nedodal, postup skončí s vysvětlením místo odhadu.

Runtime proteinový následek nepřekládá de novo z DNA sekvence. Používá toto
pořadí validovaných lokálních zdrojů:

1. normalizovaný snapshot malých indelů,
2. coding SNV snapshot,
3. ENIGMA Table 9,
4. uživatelem dodaná `p.` notace pouze tehdy, když žádný z uvedených zdrojů
   nemá validovaný následek.

Pokud lokální zdroj obsahuje kanonický následek, má tento následek přednost ve
výstupu. Uživatel může dodat také zkrácený legacy frameshift zápis, například
`p.(Cys1225fs)`. Zkrácený zápis je přijat pouze tehdy, když původní aminokyselina
a její pozice přesně souhlasí s validovaným plným následkem, například
`p.(Cys1225SerfsTer10)`. ARIANE ve výsledku vždy vrátí plný kanonický zápis.
Zápis `p.(Arg1225fs)` ani `p.(Cys1226fs)` proto pro stejnou variantu neprojde.

Příklad:

```text
BRCA1 c.303T>G
BRCA1 p.(Tyr101Ter)
variant_type = nonsense
```

Pokud uživatel zadá `p.` notaci, která není totožná ani povoleným zkráceným
frameshift ekvivalentem následku uloženého pro referenční transkript, aplikace
vrátí HTTP 422. Rozpor se neřeší výběrem jedné z hodnot.

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
| benigní | Very Strong | -8 |
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

Zdroj automatických kritérií: lokální verzovaný snapshot variantově specifických combined clinical LR odvozený z UCSC ENIGMA `BRCAmfa` tracku. Supplementary Table 7 ani její posterior probability se pro PP4/BP5 nepoužívají.

| Kritérium | Supporting | Moderate | Strong | Very Strong |
| --- | ---: | ---: | ---: | ---: |
| PP4 | LR >= 2,08 | LR >= 4,3 | LR >= 18,7 | LR >= 350 |
| BP5 | LR <= 0,48 | LR <= 0,23 | LR <= 0,05 | LR <= 0,00285 |

PP4 a BP5 se automaticky vyhodnocují z lokálního verzovaného snapshotu variantově specifických klinických LR. Manuální revize zůstává dostupná pro varianty nebo zdroje, které snapshot neobsahuje. Reviewer zadá variantově specifickou klinickou hodnotu, její škálu, citaci zdroje a souhrn zahrnutých klinických dat včetně kontroly jejich nezávislosti. Podporované škály jsou běžný LR, `log10(LR)` a ACMG evidence points. ARIANE určí sílu výhradně podle ekvivalentních prahů. Jedna publikace stačí, pokud poskytuje metodicky přijatelný variantově specifický klinický LR. Není nutné kombinovat více publikací. Sílu PP4 nelze ručně přepsat a neúplný záznam nelze aplikovat.

Automatický snapshot je uložen v souborech:

- `data/precomputed/brca_pp4_clinical_lr_snapshot.index.json`,
- `data/precomputed/brca_pp4_clinical_lr_snapshot.metadata.json`.

Snapshot je odvozen z veřejného UCSC ENIGMA `BRCAmfa` tracku verze 1.1.0. Builder `scripts/build_pp4_clinical_lr_snapshot.py` používá pouze variantově specifická data ze zdrojů uvedených v ENIGMA Appendix B, která jsou v tracku samostatně dostupná: Easton et al. 2007, PMID 17924331, Parsons et al. 2019, PMID 31131967, a Li et al. 2020, PMID 31853058. Caputo et al. 2021 je z automatického výpočtu vyřazen, protože není v použitém seznamu Appendix B. Snapshot používá přímo klinické LR. Posteriorní pravděpodobnost se nepřevádí pomocí obecného prioru.

Metadata obsahují URL a checksum zdroje, checksum indexu, verzi pravidel, počty záznamů a seznam konfliktních normalizovaných indelů. Chybějící metadata, nesprávný checksum, nesprávný počet záznamů nebo nejednoznačný alias zastaví spuštění aplikace. Pro `BRCA1 c.5266dup` a alias `c.5266dupC` obsahuje snapshot LR `6,89647 × 10^45` ze studie Li et al. 2020, což odpovídá PP4 Very Strong.

Aktuální snapshot obsahuje 4 380 jednoznačných variantových záznamů.

Zdroje uvedené v ENIGMA Appendix B jsou v aplikaci vedeny jako předem uznané metodické zdroje:

| Zdroj | PMID |
| --- | --- |
| Goldgar et al. | 15290653 |
| Thompson et al. | 12900794 |
| Easton et al. | 17924331 |
| Spurdle et al. | 25857409 |
| de la Hoya et al. | 27008870 |
| Parsons et al. | 31131967 |
| Li et al. | 31853058 |

Status zdroje se zpracovává fail-closed:

- `ENIGMA Appendix B source`: musí být vybrán jeden z uvedených PMID. Po splnění ostatních požadavků může PP4 vstoupit do amended klasifikace.
- `Other reviewed source`: vyžaduje citaci, jméno reviewera a metodické zdůvodnění kompatibility s ENIGMA PP4. Po splnění ostatních požadavků může PP4 vstoupit do amended klasifikace.
- `Unreviewed source`: hodnota a zdroj se zachovají v auditním záznamu, ale PP4 se neaplikuje a nepřidají se body.

Za primární zdroj evidence se považuje publikace nebo verzovaný dataset. Sekundární databáze, například CANVarUK, může sloužit k nalezení a zobrazení hodnoty, ale nenahrazuje citaci primárního zdroje. Pokud amended výsledek kombinuje PP1 nebo PS4 s automatickým či manuálním PP4/BP5, backend vyžaduje explicitní potvrzení nezávislosti pozorování a textové zdůvodnění. Bez nich požadavek odmítne. Samotné zaškrtnutí nenahrazuje kontrolu zdrojových kohort a klinických LR komponent.

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

Automatická Module 1 klasifikace nepřidává PS4, PM3, PP1, BS2 a BS4. PP4 a BP5 přidává pouze při přesné shodě s validovaným lokálním snapshotem klinických LR. Ostatní uvedené kódy závisejí na datech, která nelze bezpečně odvodit pouze z HGVS varianty.

Strukturovaná manuální část dále podporuje:

- PP4 z variantově specifického combined clinical LR,
- PVS1 RNA,
- BP7 RNA,
- PVS1 pro iniciační kodon,
- PS1 pro stejný splice efekt.

Manuálně doplněná kritéria vytvářejí oddělený amended working result. Původní automatická Module 1 klasifikace zůstává zachována.

## 5. Postup klasifikace

Kritéria a potřebné anotace se vyhodnocují v tomto pořadí:

1. kontrola vstupu, reference a proteinového následku,
2. souřadnice GRCh37 a GRCh38, pokud má varianta jednoznačně definovatelnou
   genomovou alelu,
3. gnomAD a BA1,
4. Table 9 pro PS3 a BS3,
5. Table 4 pro PVS1 a PM5 PTC,
6. lokální snapshot klinických LR pro PP4 nebo BP5,
7. proteinové PS1,
8. PP3 a BP4,
9. BP7,
10. BP1,
11. klasifikační kombinace.

Exonové delece a duplikace s neurčitými breakpointy nepokračují do
souřadnicového resolveru, SpliceAI ani BayesDel. Tyto služby vyžadují konkrétní
genomovou alelu a jejich použití na intervalový exonový zápis by nebylo
interpretovatelné. PVS1 se pro takovou variantu vyhodnotí přímo z ENIGMA Table 4.

BA1 ukončí klasifikaci jako Class 1.

Pokud jsou všechna aplikovaná kritéria pouze v jednom směru, používají se kombinace ENIGMA VCEP v1.2 Table 3. Samotný součet bodů v tomto případě neurčuje třídu. Tavtigian 2020 je bodový systém a pro evidenci pouze jedním směrem se v ARIANE nepoužívá.

Příklad: PVS1 Very Strong bez dalšího kritéria zůstává Class 3, protože nesplňuje kombinaci pro Likely Pathogenic.

Pokud jsou současně přítomna patogenní i benigní kritéria, nastává druhý klasifikační postup ENIGMA. V tomto případě se používá bodový systém Tavtigian 2020:

| Součet | Třída |
| ---: | --- |
| 10 a více | Class 5, Pathogenic |
| 6 až 9 | Class 4, Likely Pathogenic |
| -1 až 5 | Class 3, VUS |
| -6 až -2 | Class 2, Likely Benign |
| méně než -6 | Class 1, Benign |

Výsledek s protichůdnými směry zachovává vypočtenou ENIGMA třídu a obsahuje barevný pruh `Mixed evidence`. Pruh uvádí, že byla použita ENIGMA bodová kombinace a že je nutná expertní revize. Odkazuje přímo na verzovaný dokument [ENIGMA Specifications v1.2](https://cspec.genome.network/cspec/File/id/02537f62-66a3-4e67-8aec-cf44b326534d/data), část `Classification Methods`, druhý postup. Rozbalovací technický detail zvlášť ukazuje součet patogenních bodů, benigních bodů a celkový výsledek.

## 6. Oficiální ENIGMA tabulky

Oficiální zdroje jsou genově oddělené záznamy ClinGen CSpec pro ENIGMA BRCA1/2 VCEP v1.2.0, vydání 2025-01-09:

- BRCA1: [ClinGen CSpec GN092](https://cspec.genome.network/cspec/ui/svi/doc/GN092)
- BRCA2: [ClinGen CSpec GN097](https://cspec.genome.network/cspec/ui/svi/doc/GN097)

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

Za běhu se alias nejprve převede na kanonickou `c.` notaci. Záznam poskytne
očekávanou `p.` notaci a lokální souřadnice. Rozpor v `p.` notaci skončí chybou
422; výjimkou je pouze zkrácený frameshift zápis se shodnou původní
aminokyselinou a pozicí. Snapshot neurčuje výslednou klinickou třídu ani
automaticky nepřidává kritéria.

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

Veřejný výsledek obsahuje oddělený strukturovaný `spliceai_audit`. Hlavní klasifikace zůstává stručná a technické údaje jsou ve webovém rozhraní standardně sbalené pod položkou `Evidence details > SpliceAI`. Po rozkliknutí se zobrazí použité skóre, vybraný transkript, politika `reference_transcript`, skóre a transkript maxima přes všechny dostupné transkripty, delta pole, zdroj, GRCh38 dotaz a identifikátor cache záznamu. Stejná struktura se ukládá do auditní události dokončené klasifikace.

Zápis dynamické cache je atomický. Bez nakonfigurovaného runtime adresáře nebo
Railway volume se při lokálním vývoji používá původní `data/spliceai/`.

Intronický předpočítaný snapshot:

`data/spliceai/spliceai_brca_intronic_snv_reference_cache.json`

Aktuální snapshot je dokončen pro všech 13 800 variant souřadnicové mapy. Cache je použitelná pouze s kompletními metadaty a úspěšnou kontrolou počtu záznamů a checksumu. Rozpracovaná cache se nenačte jako platný zdroj.

### 8.3 Priorita zdrojů

Pro podporované varianty se používá lokální předpočítaná cache. Stav zdroje a důvod selhání se zaznamenávají. Chybějící skóre se nepřevádí na nulu a nemůže vytvořit BP1, BP4 ani BP7.

## 9. Souřadnice

### 9.1 Co znamená převod souřadnic

ARIANE potřebuje pro některé datové zdroje VCF-like alelu `chromosome`, `position`,
`REF` a `ALT` v GRCh37 nebo GRCh38. gnomAD v2.1.1 a MyVariant/BayesDel používají
GRCh37, zatímco gnomAD v3.1.2 a lokální SpliceAI data používají GRCh38.

Převod není prosté přičtení coding pozice ke genomové pozici. Musí respektovat:

- konkrétní verzi referenčního transkriptu,
- hranice exonů a intronů,
- orientaci genu; BRCA1 je na minus strand a BRCA2 na plus strand,
- normalizaci REF a ALT alely u inzercí, delecí a duplikací,
- rozdíly mezi GRCh37 a GRCh38.

Runtime proto souřadnice ručně nepočítá. Čte již normalizované souřadnice z
verzovaných dat nebo použije specializovaný HGVS resolver.

### 9.2 Lokální souřadnicové zdroje

Coding SNV mají GRCh37 a GRCh38 ve verzovaném coding SNV snapshotu. Malé indely
mají obě sestavy, REF a ALT ve verzovaném normalizovaném indel snapshotu.
Intronické SNV v podporovaném okně používají rozšířenou lokální souřadnicovou
mapu. Její referenční báze pocházejí z UCSC Genome Browser sequence API pro hg19
a hg38 a jsou navázány na stejné referenční transkripty jako ostatní snapshoty.

Při startu se intronická mapa a persistentní read-through cache načtou do
in-memory resolver cache. Aktuální implementace kontroluje tuto cache před
přímým lookupem coding SNV a indel snapshotu. Intronická mapa je při načítání
první a read-through cache nepřepisuje již existující klíč. Coding a indel
snapshot se použijí následně, pokud klíč nebyl nalezen v načtené resolver cache.

### 9.3 Síťové resolvery

Pro variantu, která nemá použitelné lokální souřadnice, je pořadí:

1. VariantValidator s referenčním transkriptem `NM_007294.4` nebo `NM_000059.4`,
2. Mutalyzer samostatně pro GRCh37 a GRCh38.

Výsledek má stav `ok`, pokud jsou dostupné obě sestavy, `partial`, pokud je
dostupná pouze jedna, a `failed`, pokud není dostupná žádná. Úspěšný síťový
výsledek se ukládá do read-through cache. Úplné přechodné selhání se necachuje,
aby mohl pozdější dotaz resolver zopakovat.

Nejednoznačný výsledek se nesmí převést na první nalezené ID. Exonová CNV s
neurčitými breakpointy síťové resolvery vůbec nevolá, protože nemá jednu
konkrétní VCF alelu.

Souřadnice jsou verzovaná vlastnost kombinace referenčního transkriptu,
genomového sestavení a normalizační politiky. Nemění se při každém dotazu, ale
musí se znovu validovat při změně kterékoliv z těchto částí.

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

U ostatních zdrojů se degradace ukládá do centrálního registru a zobrazuje v
klasifikačním výsledku a `/api/health`. Hlášení centrálního registru obsahuje
komponentu a konkrétní důvod. Absolutní serverové cesty se pro uživatele
zkracují na tvar `…ariane/...`.

Podrobné odpovědi síťových resolverů, HTTP těla a výjimky se zapisují do
serverového logu. Veřejný klasifikační výsledek je neopakuje. Uživatel dostane
jedno stručné a akční shrnutí pro nedostupné souřadnice, ClinVar nebo ClinGen.
Tím se zachová auditovatelná technická příčina bez zveřejnění interních cest a
dlouhých odpovědí poskytovatelů.

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

### 14.1 Regresní matice typů variant

Soubor `tests/test_variant_type_regression_matrix.py` udržuje explicitní regresní matici podporovaných typů variant. První část ověřuje odvození normalizovaného typu pro nonsense, frameshift, missense, synonymous, intronickou SNV, canonical splice-site variantu, in-frame indel, exonovou deleci a duplikaci, iniciační kodon a obě UTR oblasti.

Druhá část pro reprezentativní varianty kontroluje celý deterministický průchod Module 1 od typu varianty přes lokální tabulky a snapshoty po aplikovaná kritéria a výslednou třídu. Každý řádek obsahuje očekávaná i zakázaná kritéria. Tím se kontroluje nejen přítomnost správné větve, ale také nepřípustné přičtení kritéria z jiné větve, například PP3 u nonsense nebo frameshift varianty.

Synonymous varianty mají v matici obě ENIGMA větve. Varianta uvnitř funkční domény s potvrzeným SpliceAI nejvýše 0,1 vyžaduje BP4 a BP7 a nesmí dostat BP1. Varianta mimo funkční domény se stejným limitem SpliceAI vyžaduje BP1 Strong a nesmí současně dostat BP4 ani BP7. Pozitivním regresním případem pro BP1 je `BRCA1 c.306A>C p.(Ala102=)`.

RNA-dependent canonical splice-site a initiation-codon větve samostatně ověřují, že aplikace nevytvoří automatické PVS1 a místo něj vrátí doporučení ke strukturované manuální revizi. Matice nepoužívá živé externí služby. Používá verzované lokální ENIGMA tabulky a snapshoty, takže musí být deterministická.

### 14.2 Interakce funkční, RNA a predikční evidence

ARIANE nepoužívá plošné potlačení podle příznaku `has_functional_evidence`. Interakce se vyhodnocují podle mechanismu a hierarchie důkazů z ENIGMA v1.2 Figure 1A, Figure 1B, Figure 1C a Appendix E.

Přijaté PVS1 (RNA) nahrazuje slabší bioinformatické kódy pro stejný experimentálně potvrzený splice důsledek. Přijaté BP7 Strong (RNA) nahrazuje BP7 Supporting, ale podle Figure 1B obecně zachovává ostatní použitelné bioinformatické kódy. PS3 nebo BS3 bez PVS1 automaticky nepotlačuje PP3, BP4, BP7 ani BP1, protože Figure 1C výslovně požaduje zachování relevantních bioinformatických kódů.

Každé nahrazení, zachovaná potenciální interakce nebo konflikt se vrací ve strukturovaném poli `evidence_interactions`. Ve webovém rozhraní je zobrazeno v rozbalovací části `Evidence interaction warnings`. Přesná matice je v `docs/evidence_interaction_matrix.md`.

Související dokumenty:

- `docs/enigma_source_data_audit.md`,
- `docs/evidence_interaction_matrix.md`,
- `docs/manual_evidence_review.md`,
- `docs/splice_ps1_reference_set.md`,
- `docs/vus_explanation_and_golden_cases.md`.
