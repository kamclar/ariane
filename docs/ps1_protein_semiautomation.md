# Proteinové PS1 v ARIANE

## Rozhodovací pravidlo

Proteinový PS1 registr spojuje P/LP missense reference z ENIGMA ST7 v1.2 s
aktuálními assertion ENIGMA BRCA1/2 VCEP v1.2 publikovanými v ClinGen Evidence
Repository. Aktuální registr obsahuje 85 referencí: 63 `eligible` a 22
`excluded` podle známé RNA a splice evidence. Z nich 26 pochází z aktuálních
ERepo v1.2 assertion a 60 ze ST7; jedna reference je v obou zdrojích. Nejde o
bezpodmínečný PS1 allowlist. SpliceAI skóre není v registru uloženo a kontroluje
se při každém použití reference.

Automatické proteinové PS1 vyžaduje:

1. stejnou normalizovanou missense substituci a jinou nukleotidovou změnu;
2. P/LP klasifikaci reference v ENIGMA ST7 v1.2, aktuální ENIGMA/ClinGen VCEP
   assertion nebo úplné lokální reklasifikaci podle uvedené verze pravidel;
3. SpliceAI nejvýše 0,1 u reference i hodnocené varianty;
4. žádný známý škodlivý splice efekt po kontrole uvedených verzovaných zdrojů;
5. stav `eligible` v proteinovém PS1 registru.

SpliceAI pro referenci i hodnocenou variantu se získá při klasifikaci ze stejné
nakonfigurované služby. Chybějící skóre se nepovažuje za nulu a PS1 se
automaticky nepřidělí.

Implementační politika ARIANE používá oficiální P/LP klasifikaci v ENIGMA ST7
v1.2 jako přijatý klasifikační základ reference. Toto rozhodnutí bylo potvrzeno
odbornou metodickou konzultací 7. září 2026. ST7 klasifikace vznikly v
historickém ENIGMA multifaktoriálním likelihood referenčním souboru. Nešlo o
aplikaci proteinového PS1, proto je u těchto záznamů závislost na PS1 vedena
jako `false`.

P reference dává PS1 Strong. LP reference dává PS1 Moderate. Síly z více
referencí se nesčítají.

```text
ST7 kandidát
          |
          v
stejný missense následek + jiná c. změna
          |
          v
přijatá P/LP reference ze ST7 nebo VCEP?
     |                         |
    ne                        ano
     |                         |
bez PS1                  splice podmínky splněny?
                          |             |
                         ne            ano
                          |              |
                    bez PS1 / revize  automatické PS1
                                      P: Strong
                                      LP: Moderate
```

## Datasety

- `st7_reference_set.json`: úplný oficiální zdroj kandidátů;
- `enigma_st2_splice_evidence.json`: úplných 220 řádků ENIGMA ST2 a všech 383
  navázaných referencí ST3 pro kontrolu a doložení známých RNA výsledků;
- `enigma_table9.json`: funkční a publikovaná splice evidence;
- `enigma_erepo_vcep_registry.json`: úplný lokální snapshot 180 publikovaných
  assertion panelu ENIGMA BRCA1/2 VCEP, včetně přesné verze assertion metody;
- `enigma_erepo_vcep_registry.metadata.json`: checksum, datum získání, verze API
  a počty záznamů podle verze a stavu;
- `ps1_protein_reference_registry.json`: 85 P/LP missense referencí z aktuálního
  ERepo v1.2 a ST7 s explicitním stavem a auditními podklady;
- splice PS1 nemá aktivní referenční registr; vyžaduje samostatnou strukturovanou manuální revizi.

## Co patří do registru

Registr přijímá pouze missense P/LP reference v kanonickém ENIGMA transkriptu,
které mají dohledatelný klasifikační původ. Přípustné klasifikační základy jsou:

1. P/LP kandidáti z oficiální ENIGMA Supplementary Table 7 v1.2;
2. aktuální verzovaná ENIGMA BRCA1/2 VCEP v1.2 assertion z checksumovaného
   lokálního snapshotu ClinGen Evidence Repository;
3. lokální úplná reklasifikace podle uvedené verze ENIGMA VCEP pravidel.

Lokální reklasifikace musí být označena
`locally_recurated_under_enigma_vcep`. Nesmí se vydávat za oficiální expert
panel assertion. Musí obsahovat identifikátor posuzovatele, datum, verzi
pravidel a identifikátor úplného evidenčního záznamu.

Nové oficiální ERepo reference se získají reprodukovatelným builderem. Lokálně
reklasifikované reference mimo ERepo a ST7 se zapisují do
`backend/data/ps1_protein_reference_extensions.json`. Generátor všechny zdroje
sloučí do jediného runtime registru.

Při nálezu reference ARIANE předvyplní referenční c. a p. notaci, přijatou
klasifikaci a zdroj, shodu proteinového následku, rozdílnou nukleotidovou změnu,
dostupná SpliceAI skóre a výsledky kontroly definovaných RNA/splice zdrojů.
Aktuální ERepo v1.2 assertion má přednost před historickým zdrojem ST7.
ClinVar a živé ERepo slouží jen jako externí porovnání. Uživatel potvrzuje jen
podmínky, které nebylo možné doložit automaticky. Neúplná revize nepřidá body.

ClinVar aggregate, počet hvězdiček, CANVarUK, BRCA Exchange, jednotlivá
publikace nebo výpočetní predikce nestačí k vytvoření `eligible` reference.
Tříhvězdičková historická ENIGMA assertion se zobrazí s upozorněním, ale bez
odpovídající aktuální v1.2 assertion v lokálním ERepo registru se automaticky
nepoužije.

Ve strukturované manuální revizi lze zadat pouze c. HGVS referenční varianty.
Backend z referenčního transkriptu odvodí a ověří p. následek, porovná jej s
hodnocenou variantou a získá SpliceAI pro obě varianty. Přesnou referenci ověří
také v lokálním ERepo registru a pro porovnání v ClinVar. Hvězdičky ClinVar popisují pouze review status
a žádný jejich počet nekvalifikuje referenci pro PS1. Běžný ClinVar aggregate
závěr nepředvyplňuje klasifikaci, ověření, klasifikační zdroj ani evidenční
reference. Použít lze pouze samostatně identifikovanou assertion příslušného
ENIGMA/ClinGen VCEP v lokálním ERepo snapshotu nebo úplnou lokální reklasifikaci.

## Odkud se berou jednotlivá pole

| Informace v registru | Zdroj |
|---|---|
| gen, c. a p. notace, P/LP třída a původní klasifikační zdroj | ST7; u rozšíření oficiální VCEP assertion nebo lokální evidenční záznam |
| referenční transkript | BRCA1 `NM_007294.4`, BRCA2 `NM_000059.4` podle ENIGMA v1.2 |
| normalizovaná missense substituce | kanonická p. notace ověřená proti referenčnímu transkriptu |
| podklad proteinového mechanismu | PS3 funkční evidence z Table 9, nebo u patogenní missense reference doložená absence predikovaného a potvrzeného splice efektu |
| SpliceAI reference | výpočet na požádání stejnou profilově připnutou službou jako u hodnocené varianty |
| známá RNA/splice evidence | úplná ENIGMA Table 9 v1.2 a úplná Supplementary Table 2 v1.2 |
| stav `eligible`, `excluded`, `review_required` | ST7 je `eligible`, pokud definované zdroje nevylučují proteinovou větev; známý škodlivý splice efekt vede k `excluded`; konflikt nebo neúplný podklad vede k `review_required`; predikční podmínka se ověřuje za běhu |
| případná PS1 závislost klasifikace reference | u ST7 historický multifaktoriální klasifikační postup; u nové assertion nebo lokální reklasifikace úplný evidenční záznam |
| provenance a checksumy | generátor registru ze všech použitých verzovaných vstupních souborů |

Registr nyní obsahuje 40 referencí se stavem `eligible`. Z nich má 35 v Table 9
PS3 Strong funkční evidenci. U zbývajících
pěti je proteinový mechanismus zaznamenán jako patogenní missense reference bez
predikovaného nebo potvrzeného splice efektu. Přímý proteinový funkční test tedy
není povinný pro každý záznam, ale použitý mechanistický podklad musí být v
registru explicitní.

Stav `eligible` vyžaduje uzavřenou kontrolu závislosti klasifikace reference na
PS1. U ST7 je uvedeno `false`, protože zaznamenaná IARC klasifikace pochází z
multifaktoriálního likelihood referenčního souboru, nikoliv z aplikace
proteinového PS1. U novějších assertions musí evidenční záznam uvést, zda bylo
PS1 použito. Pokud ano, musí uvést použité reference a validátor vyloučí přímou
i delší známou kruhovou závislost.

Runtime SpliceAI záznam ukládá profil, referenční genom, transkript, vstupní
variantu, parametry výpočtu a zdroj. Změna modelu nebo anotace vytváří nový
profil a starší runtime záznam se nepoužije.

Pro hodnocenou i referenční variantu se používá výsledek stejného
ENIGMA kompatibilního SpliceAI zdroje. Hodnota `spliceai_prediction` uvedená v
Table 9 popisuje kontext, ve kterém ENIGMA posoudila funkční evidenci PS3/BS3.
Nepřepisuje aktuální predikční výsledek a při jeho nedostupnosti neslouží jako
fallback pro proteinové PS1. Rozdíl se zaznamená do auditu. Pokud hodnoty leží
v různých ENIGMA pásmech, je nutná odborná kontrola jejich provenance.

`none_identified` znamená pouze, že v uvedených verzích definovaných zdrojů
nebyl nalezen odpovídající škodlivý splice záznam. Nejde o tvrzení, že žádná
RNA evidence neexistuje.

Registr se ověří příkazem:

```bash
python scripts/validate_ps1_protein_registry.py
```

Registr se reprodukovatelně vytvoří příkazem:

```bash
python scripts/build_ps1_protein_reference_registry.py
```

Před sestavením PS1 registru se aktualizuje ERepo snapshot:

```bash
python scripts/build_enigma_erepo_vcep_registry.py
```

Builder stáhne úplný ERepo export, vybere přesně panel `ENIGMA BRCA1 and BRCA2
VCEP` a u každého záznamu načte detail assertion metody. Pouze verze `1.2.0`
má stav `current_vcep_assertion`. Starší verze a záznamy bez uvedené verze jsou
uloženy pro audit a upozornění, ne pro automatickou způsobilost PS1.

Po úmyslné kurátorské změně schváleného záznamu lze jeho checksum přepočítat
výhradně explicitním příkazem:

```bash
python scripts/validate_ps1_protein_registry.py --write-checksums
```

Přepočítání checksumu samo nenahrazuje odborné schválení obsahu.

## Kruhová závislost

Pokud je známo, že klasifikace reference použila PS1, musí být uložena použitá
reference. Bez ní záznam nelze schválit. Validátor odmítá přímou závislost na
sobě a známé cykly. Způsobilá nová assertion musí mít stav závislosti `true`
nebo `false`; neznámý stav zůstává k revizi. U ST7 je `false` odvozeno z
multifaktoriálního klasifikačního postupu zaznamenaného v oficiálním Appendixu.
