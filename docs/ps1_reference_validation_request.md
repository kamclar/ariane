# Potvrzení podmínek proteinového PS1

## Stav rozhodnutí

Vyřešeno 7. září 2026. Odborná metodická konzultace potvrdila, že P/LP
zařazení v ENIGMA ST7 v1.2 lze přijmout jako splnění klasifikačního požadavku
reference pro proteinové PS1. ST7 slouží jako ochrana proti převzetí běžné
historické ClinVar klasifikace bez ENIGMA podkladu.

Implementace přijímá ST7 pouze jako klasifikační základ reference. Nadále
vyžaduje stejnou normalizovanou missense substituci z jiné nukleotidové změny,
proteinový mechanismus, SpliceAI nejvýše 0,1 pro obě varianty a kontrolu známé
RNA/splice evidence. Po přidání aktuálního ERepo v1.2 snapshotu registr obsahuje
85 referencí: 63 záznamů `eligible` a 22 záznamů `excluded`. Původních 60 ST7
referencí zůstává dohledatelných přes `source_memberships`.

ST7 uvádí IARC klasifikace z historického multifaktoriálního likelihood
referenčního souboru. Tyto klasifikace nevznikly použitím proteinového PS1,
proto je jejich PS1 závislost zaznamenána jako `false`. Novější VCEP assertion
se při zobrazení nadále upřednostní a její použitá kritéria se uchovají pro
audit.

Následující část zachovává původní dotaz pro audit rozhodnutí.

## Účel

ARIANE měla 60 P/LP missense referencí převzatých z ENIGMA Supplementary
Table 7 v1.2. Před metodickým rozhodnutím bylo 40 vedeno jako
`review_required` a 20 jako `excluded` kvůli známému splice efektu.

Ze 40 kandidátů k revizi je 31 pro BRCA1 a 9 pro BRCA2. U 35 je v Table 9
zaznamenána PS3 Strong funkční evidence. U pěti je proteinová větev podložena
missense mechanismem bez nalezeného predikovaného nebo potvrzeného splice
efektu v definovaných zdrojích.

## Navrhovaná podmínka pro stav `eligible`

Reference by mohla být použita pro automatické proteinové PS1 pouze při
současném splnění všech následujících podmínek:

1. Jde o missense variantu v referenčním transkriptu daného genu.
2. Reference je Pathogenic nebo Likely Pathogenic podle příslušných VCEP
   specifications.
3. Hodnocená a referenční varianta vytvářejí stejnou normalizovanou
   aminokyselinovou substituci, ale vznikají jinou nukleotidovou změnou.
4. SpliceAI je nejvýše 0,1 u obou variant.
5. U žádné z obou variant není v definovaných ENIGMA zdrojích potvrzen škodlivý
   splice efekt.
6. Je ověřeno, zda klasifikace reference použila PS1. Neznámý stav nepřidělí
   body.
7. Pokud klasifikace reference PS1 použila, jsou uvedeny její PS1 reference a
   je vyloučena přímá i delší kruhová závislost.
8. Každý podklad má stabilní identifikátor, verzi pravidel, datum kontroly a
   auditní odkaz.

Při splnění podmínek by Pathogenic reference vedla k PS1 Strong a Likely
Pathogenic reference k PS1 Moderate.

## Původní otázky k metodickému potvrzení

1. Stačí P/LP zařazení varianty v ST7 jako doklad, že byla klasifikována podle
   VCEP specifications, nebo je potřeba samostatná ENIGMA/ClinGen VCEP
   assertion?
2. Lze proteinový mechanismus pro tento účel potvrdit absencí predikovaného a
   potvrzeného splice efektu při SpliceAI nejvýše 0,1, nebo je vždy nutná také
   přímá proteinová funkční evidence?
3. Může být použita reference, jejíž klasifikace sama obsahovala PS1, pokud jsou
   známé všechny závislosti a je prokázáno, že nevzniká kruh? Nebo mají být
   takové reference pro automatické PS1 zcela vyloučeny?
4. Který oficiální zdroj je nejvhodnější pro ověření klasifikace reference,
   mechanismu a použitých kritérií? ClinGen Evidence Repository, samostatná
   ENIGMA expert-panel assertion v ClinVaru, nebo jiný ENIGMA soubor?
5. Může ENIGMA poskytnout seznam schválených proteinových PS1 referencí včetně
   třídy, verze pravidel, mechanismu a případných PS1 závislostí?
6. Je správné odvodit PS1 Strong z Pathogenic reference a PS1 Moderate z Likely
   Pathogenic reference po splnění všech výše uvedených podmínek?

## Příklad k ověření

Pro hodnocenou variantu `BRCA1 NM_007294.4:c.5217T>A p.(Asp1739Glu)` je
potenciální referencí jiná nukleotidová změna se stejným proteinovým následkem,
`BRCA1 NM_007294.4:c.5217T>G p.(Asp1739Glu)`. Tato reference není v současném
ARIANE ST7 registru. Potřebujeme ověřit, zda má platnou ENIGMA VCEP P/LP
klasifikaci a zda splňuje všechny podmínky proteinového PS1.

## Text k odeslání Janě bez diakritiky

Ahoj Jano,

chteli bychom potvrdit postup pro automaticke proteinove PS1. ARIANE ma ze ST7
40 P/LP missense kandidatu, ktere po kontrole znamych RNA a splice dat mohou
patrit do proteinove vetve. Zatim ale zadny z nich automaticky nepovolujeme.

Navrhujeme povolit konkretni referenci pouze tehdy, kdyz dolozime P/LP
klasifikaci podle prislusnych VCEP specifications, stejny normalizovany missense
nasledek z jine nukleotidove zmeny, SpliceAI nejvyse 0,1 u reference i hodnocene
varianty a zadny potvrzeny skodlivy splice efekt u obou variant. Zaroven bychom
vzdy overili, zda klasifikace reference pouzila PS1. Pokud ano, zaznamenali
bychom pouzitou referenci a vyloucili primou i delsi kruhovou zavislost. Neznamy
stav by zadne body nepridelil.

Muzeme se prosim zeptat:

1. Staci zarazeni P/LP varianty v ST7 jako doklad klasifikace podle VCEP
   specifications, nebo potrebujeme samostatnou ENIGMA/ClinGen VCEP assertion?
2. Staci pro potvrzeni proteinoveho mechanismu absence predikovaneho a
   potvrzeneho splice efektu pri SpliceAI nejvyse 0,1, nebo je vzdy nutna i prima
   proteinova funkcni evidence?
3. Muzeme pouzit referenci, jeji klasifikace sama obsahovala PS1, pokud zname
   vsechny zavislosti a prokazeme, ze nevznika kruh? Nebo mame takove reference
   zcela vyloucit?
4. Ktery oficialni zdroj mame pouzit pro overeni klasifikace, mechanismu a
   pouzitych kriterii reference?
5. Je spravne po splneni vsech podminek priradit PS1 Strong podle Pathogenic
   reference a PS1 Moderate podle Likely Pathogenic reference?

Jako konkretni priklad resime BRCA1 NM_007294.4:c.5217T>A
p.(Asp1739Glu). Potencialni reference BRCA1 NM_007294.4:c.5217T>G vede ke
stejnemu proteinovemu nasledku, ale neni v nasem soucasnem ST7 registru. Ma tato
varianta platnou ENIGMA VCEP P/LP klasifikaci a lze ji po uvedenych kontrolach
pouzit jako referenci pro proteinove PS1?

Pokud existuje oficialni seznam vhodnych PS1 referenci vcetne pouzitych
kriterii nebo zavislosti, radi bychom jej pouzili jako verzovany podklad.

Dekuji.

## Přiložený pracovní přehled

Soubor `ps1_reference_review_candidates.tsv` se generuje z aktuálního registru.
Po metodickém rozhodnutí neobsahuje žádný ST7 záznam čekající pouze na potvrzení
klasifikačního základu. Původní seznam 40 kandidátů zůstává dohledatelný v
historii projektu.
