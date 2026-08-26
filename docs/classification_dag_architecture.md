# DAG architektura klasifikace ARIANE

## Účel

Migrace odděluje normalizaci vstupu, získání evidence, rozhodování o kritériích,
řešení interakcí evidence, klasifikaci a prezentaci. Přesun výpočtu do grafu
nesmí sám o sobě změnit klinický výsledek.

Produkční verze `4.0.0-gene-policy-provider-dag` již starý sekvenční klasifikátor nevolá.
Vytváří `NormalizedVariant`, `EvidenceBundle`, samostatné výsledky rodin kritérií,
mechanism-aware interakce, výslednou ENIGMA klasifikaci a manuální review výstupy.
Každá rodina má vlastní uzel a společný slovník výsledku vzniká až v poslední
projekci do veřejného API.

Pravidlové uzly jsou rozdělené podle odpovědnosti v
`backend/classification_dag/nodes/`:

| Modul | Odpovědnost |
| --- | --- |
| `context.py` | SpliceAI kontext a porovnání provenance s Table 9 |
| `population.py` | BA1, BS1, PM2 a populační větev exonových CNV |
| `functional.py` | Kalibrovaná funkční evidence PS3 a BS3 z Table 9 |
| `pvs1.py` | PVS1, PVS1 RNA a PM5 PTC |
| `clinical.py` | Klinické LR PP4/BP5 a proteinové PS1 |
| `bioinformatic.py` | Figure 1A, PP3, BP4, BP1 a BP7 |
| `policy.py` | Interakce evidence a výsledná ENIGMA kombinace |
| `review.py` | Vytvoření požadavků na manuální revizi |
| `support.py` | Neměnné společné typy a převod rozhodnutí |

Produkční `runtime.py` importuje tematické moduly přes `nodes/__init__.py`.
Původní `native.py` je pouze kompatibilní re-export a neobsahuje rozhodovací
logiku. Je dočasný a nesmí se stát produkční závislostí. Po přijetí DAG migrace
se odstraní společně se starým `backend/modules/classifier.py` a paritními testy.

Klasifikační evidence se již nezískává v `main.py`. Samostatné provider uzly
získávají souřadnice, SpliceAI, BayesDel, gnomAD, Table 9, klinické LR, exonové
CNV, důležité aminokyselinové pozice a podklady proteinového PS1. Každý zdroj
vrací `EvidenceItem` s explicitním stavem, důvodem a provenance. Sestavovací uzel
vytvoří jediný neměnný `EvidenceBundle`, ze kterého čtou pravidlové uzly.

## Aplikační orchestrace

FastAPI modul `backend/main.py` je kompoziční kořen a transportní vrstva. Neřeší
normalizaci, výběr providerů, spuštění klasifikačního DAGu, diagnostiku zdrojů ani
sestavení veřejného klasifikačního modelu.

Tyto odpovědnosti jsou rozdělené v `backend/services/`:

| Služba | Odpovědnost |
| --- | --- |
| `evidence_orchestration.py` | `ClassificationCommand`, normalizace, typ varianty, sestavení `ClassificationRequest`, paralelní provider DAG a externí porovnání, diagnostika dostupnosti |
| `classification_presentation.py` | Převod strukturovaného výsledku a artefaktů na stabilní `ClassificationResult`, bez změny kritérií nebo třídy |
| `variant_classification_service.py` | Funkce `execute_variant_classification()`, jediný aplikační use case pro klasifikaci jedné varianty, který spojuje orchestraci a prezentaci |

Orchestrátor přijímá volitelné `ProviderDependencies` a
`ExternalEvidenceDependencies`. Produkce používá adaptéry registrované v DAGu,
testy mohou dodat deterministické implementace. Výběr genu a VCEP politiky stále
probíhá přes verzovaný gene policy manifest. Přidání genu proto nemění transportní
ani orchestration vrstvu.

ClinVar a ClinGen ERepo nejsou vstupem klasifikace. Běží paralelně s DAGem a
připojují se až jako externí porovnání vypočteného výsledku.

## Hranice vrstev

```text
zadaná notace
      |
      v
normalizace a kontrola reference
      |
      v
NormalizedVariant
      |
      v
DAG získání evidence
      |
      v
EvidenceItem
      |
      v
DAG rozhodování o kritériích
      |
      v
CriterionDecision
      |
      v
interakce evidence podle mechanismu
      |
      v
klasifikační politika ENIGMA
      |
      v
prezentace a externí porovnání
```

Normalizace vstupu není klasifikační uzel. Porovnání s ClinVar a ClinGen ERepo
je následná anotace a nesmí ovlivnit vypočtenou třídu.

## Kontrakt uzlu

Každý uzel deklaruje neměnný identifikátor a verzi, požadované vstupní klíče,
poskytované výstupní klíče a vyhodnocovací funkci. Executor kontroluje, že:

1. každý požadovaný vstup má právě jednoho providera nebo je vstupem grafu;
2. každý výstup má právě jednoho providera;
3. graf neobsahuje cyklus;
4. úspěšný uzel vrátí přesně deklarované výstupy;
5. neúspěšný stav nepublikuje hodnoty;
6. výjimka nebo porušení kontraktu zastaví výpočet a označí chybný uzel.

Stavy uzlu jsou `succeeded`, `not_applicable`, `unavailable`, `ambiguous`,
`review_required`, `failed` a `skipped`. Musí zůstat významově odlišné.
Nedostupná data se zejména nesmí převést na nesplněné kritérium.

## Provozní režim

Proměnná `ARIANE_CLASSIFIER_ENGINE` přijímá pouze `dag`, který je zároveň výchozí
hodnotou. Hodnoty `legacy`, `shadow` ani tichý fallback nejsou povoleny.

Aktivní režim vrací `/api/health`. DAG zapisuje do logu strukturovanou stopu všech
uzlů. Porovnání původní a DAG cesty patří pouze do testovacího běhu, aby
vývojová kontrola nezvyšovala dobu odezvy a nemohla zasáhnout uživatele.

Asynchronní executor spouští nezávislé provider uzly paralelně. Blokující
knihovny a HTTP klienti běží mimo event-loop vlákno. Porucha uzlu nebo porušení
jeho kontraktu zastaví klasifikaci. Očekávaná nedostupnost externího skóre je
naopak platný výsledek provideru s `EvidenceStatus.UNAVAILABLE`; pravidla pak
nedostupnou hodnotu nesmějí zaměnit za nulové skóre nebo nesplněné kritérium.
Auditní záznam obsahuje pro každý `EvidenceItem` stav, identifikaci zdroje,
verzi, checksum a důvod. Podrobná zdrojová provenance zůstává oddělená od
stručného klinického výsledku.

## Invarianty migrace

1. Pravidla přímo neotevírají JSON soubory ani nevolají API. Přijímají typovanou
   evidenci od provider uzlů.
2. Oficiální tabulky ENIGMA jsou verzované zdroje evidence, nikoli fallbacky.
3. Žádné kritérium se nepřiřazuje podmínkou pro jednu konkrétní variantu.
4. Každé použití, vyloučení, nedostupnost nebo požadavek na revizi zaznamená
   verzi pravidla, verzi a checksum zdroje a důvod rozhodnutí.
5. Klasifikační uzel pracuje pouze se zachovanými rozhodnutími o kritériích a
   neprovádí externí lookup.
6. Textové vysvětlení vzniká ze strukturovaných rozhodnutí a neobsahuje vlastní
   paralelní klasifikační logiku.
7. Nový gen se přidává do checksumovaného
   `backend/data/gene_policy_manifest.json`. Záznam určuje transkript, protein,
   VCEP policy, prahy, funkční domény, použitelná pravidla a požadovaná data.
   Providery a regresní data musí pokrýt deklarované `required_rule_data`.
8. Pole `implementation_profile` váže manifest na konkrétní implementaci DAGu.
   Neznámý profil se odmítne před získáním evidence. Žádný nový gen nedědí BRCA
   rozhodovací graf pouze proto, že používá podobně pojmenované ACMG kódy.
9. Genové prefixy ve vstupu, zdrojové CSpec URL, HGVS startup testy, PVS1 obrázky
   a popisy funkčních domén se načítají z manifestu. Python neobsahuje větvení
   `if gene == "BRCA1"` nebo `if gene == "BRCA2"`.

## Implementované uzly

1. vstupní kontrakt a normalizovaná identita varianty;
2. provider souřadnic GRCh37 a GRCh38;
3. provider SpliceAI včetně skóre kandidátních PS1 referencí;
4. provider BayesDel a informační anotace AlphaMissense;
5. provider gnomAD;
6. provider Table 9 s verzí a SHA-256;
7. provider klinických LR, exonových CNV a důležitých reziduí;
8. provider proteinového PS1;
9. sestavení a validace `EvidenceBundle`;
10. populační BA1, BS1 a PM2;
11. Table 9 PS3 a BS3;
12. PVS1, PVS1 RNA a PM5 PTC;
13. klinické LR PP4 a BP5;
14. proteinové PS1;
15. Figure 1A PP3, BP4, BP7 a BP1;
16. interakce a deduplikace evidence;
17. ENIGMA Table 3 nebo bodová metoda pro mixed evidence;
18. manuální review triage;
19. validace `VariantAssertion` a veřejná projekce.

Uzel `policy.enigma_combination` nerozhoduje o výjimce jediného Strong
benigního kritéria podle názvu kódu. Čte explicitní vlastnost
`single_strong_likely_benign_eligible`, kterou smí vytvořit pouze příslušný
pravidlový uzel z auditovatelné provenance. BP1 ji získá po úplném průchodu
Figure 1A. BP5 ji získá z více LR příspěvků a klinických typů v kontrolovaném
datasetu. Manuální BS4 ji získá pouze z více strukturovaných LR komponent s
unikátními skupinami nezávislosti. Chybějící vlastnost se vyhodnotí jako
nezpůsobilá.

Samostatný graf `ariane.vcep.manual-evidence` zpracovává odborně doplněná
kritéria. Odděluje vstupní validaci, výpočet síly manuálních kritérií,
mechanism-aware interakce a výslednou ENIGMA kombinaci. Chyba kteréhokoli uzlu
zastaví výpočet a nevydá review-adjusted klasifikaci.

BA1 zůstává stand-alone benigní klasifikací také po doplnění manuální
evidence. Ostatní kritéria zůstávají ve výsledku pro audit, ale nepřepínají
klasifikaci do bodového postupu pro mixed evidence.

Starý `backend/modules/classifier.py` není produkční závislost ani dostupný runtime
režim. Po samostatné akceptaci DAG a přesunu zbývajících testů se odstraní tento
soubor, paritní testy a dočasný kompatibilní `backend/classification_dag/native.py`.

Synchronní graf `4.0.0-gene-policy`, který přijímá předem sestavené
`ClassificationInputs`, zůstává pouze pro izolované testy pravidel a porovnání se
starým testovacím oraclem. HTTP klasifikace vždy používá asynchronní provider graf
`4.0.0-gene-policy-provider-dag`.

## Vazba na externí modely

Interní model zůstává malý a stabilní, ale používá stejné oddělení evidence,
evidence lines a výsledných tvrzení jako ClinGen a GA4GH. Případný export do
VA-Spec bude samostatný adaptér, aby runtime nezávisel na vyvíjejícím se
externím schématu.

Použité architektonické reference:

- ClinGen Interpretation Model:
  <https://dataexchange.clinicalgenome.org/interpretation/>
- GA4GH VA-Spec Evidence Line:
  <https://va-spec.ga4gh.org/en/stable/core-information-model/entities/information-entities/evidence-line.html>
- biocommons hgvs:
  <https://pmc.ncbi.nlm.nih.gov/articles/PMC6282708/>
- VariantValidator:
  <https://github.com/openvar/variantValidator>

## Společná vrstva pravidel a uživatelský audit

Strojově čitelný katalog `backend/data/enigma_rule_catalog.json` je společným
zdrojem pro DAG, veřejné API a zobrazení pravidel. Obsahuje identifikaci ENIGMA
GN092 v1.2.0, oficiální odkazy a SHA256 zdrojových dokumentů. Uzlové a
geometrické definice dalších diagramů jsou v
`backend/data/enigma_rule_diagrams.json`. Podmínky nejsou spustitelné řetězce.
Každý `condition_id` odkazuje na konkrétní verzovanou implementaci v Pythonu.

Figure 1A je v katalogu zachycena ve všech třech větvích: missense/in-frame,
synonymous a intronic. Kritéria PP3, BP4, BP7 a BP1 již při vlastním vyhodnocení
vytvářejí strukturovanou `decision_path`. Ta obsahuje navštívené uzly, výsledek
každého testu, pozorovanou hodnotu, výsledný uzel a oficiální zdroj. Nejde o
pozdější rekonstrukci z textového důvodu.

Veřejné endpointy jsou:

- `/api/rules` pro verze, zdroje a seznam dostupných pravidel;
- `/api/rules/trees/{tree_id}` pro každý validovaný oficiální nebo odvozený
  diagram;
- `/api/rules/tables/table9` pro filtrovaný a stránkovaný veřejný pohled na
  runtime Table 9;
- `/api/rules/tables/{table_id}` pro stránkovaný a filtrovaný pohled na všech
  42 Specification, Appendix a Supplementary Tables.

API nevydává lokální cesty ani celé zdrojové sešity. Table 9 ve specializovaném
endpointu používá stejný validovaný JSON jako klasifikátor. Obecný prohlížeč
používá oddělený, reprodukovatelný prezentační snapshot všech tabulek. Hlavní
výsledek zůstává stručný;
rozhodovací cesta je rozbalovací detail u aplikovaného kritéria a celý strom je
na samostatné stránce ENIGMA rules.

Diagram se vykresluje jako SVG se směrovanými hranami. U oficiálních ENIGMA
Figures používá ruční layout podle předlohy; odvozené diagramy mají v API i UI
provenienci `ariane_derived`. `decision_path`
zvýrazňuje navštívené uzly a hrany a do uzlů doplňuje skutečné podklady
rozhodnutí, například hodnotu SpliceAI, BayesDel nebo určenou funkční doménu.
Původní panely Figure 1A, 1B a 1C a Appendix Figures 1 až 6 a 9 jsou navíc
vyjmuty z ověřených zdrojových dokumentů do samostatných statických obrazových
souborů. Gene maps Figures 1/2 současně nesou podklady používané ve Figures 7/8.
Odkaz na obrázek proto otevírá originální panel, zatímco oddělený odkaz otevírá
úplný oficiální dokument.
