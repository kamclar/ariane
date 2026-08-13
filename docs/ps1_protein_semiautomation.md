# Proteinové PS1 v ARIANE

## Rozhodovací pravidlo

ST7 je oficiální ENIGMA referenční dataset. Jeho P/LP missense varianty jsou
zařazeny do proteinového PS1 registru. Aktuální registr obsahuje 60 referencí:
35 `eligible` a 25 `excluded`.

Automatické proteinové PS1 vyžaduje:

1. stejnou normalizovanou missense substituci a jinou nukleotidovou změnu;
2. P/LP klasifikaci reference ověřenou podle ENIGMA/ClinGen VCEP;
3. SpliceAI nejvýše 0,1 u reference i hodnocené varianty;
4. žádný známý škodlivý splice efekt po kontrole uvedených verzovaných zdrojů;
5. stav `eligible` v proteinovém PS1 registru.

Implementační politika ARIANE přijímá oficiální P/LP klasifikaci v ENIGMA ST7
v1.2 jako klasifikační základ reference. ST7 záznam se však stane `eligible`
teprve po samostatné kontrole proteinového typu a splice podmínek.

P reference dává PS1 Strong. LP reference dává PS1 Moderate. Síly z více
referencí se nesčítají.

```text
ST7 a proteinový PS1 registr
          |
          v
stejný missense následek + jiná c. změna
          |
          v
stav reference v registru?
     |                  |
review/excluded       eligible
     |                  |
revize nebo důvod   SpliceAI VUA <= 0,1
                    a bez známého splice efektu?
                         |              |
                        ne             ano
                         |              |
                  bez PS1 / revize   automatické PS1
                                     P: Strong
                                     LP: Moderate
```

## Datasety

- `st7_reference_set.json`: úplný oficiální zdroj kandidátů;
- `enigma_st2_splice_evidence.json`: úplných 220 řádků ENIGMA ST2 pro
  kontrolu známých RNA výsledků;
- `enigma_table9.json`: funkční a publikovaná splice evidence;
- `ps1_protein_reference_registry.json`: všech 60 P/LP missense referencí ST7
  s explicitním stavem a auditními podklady;
- `splice_ps1_reference_set.json`: samostatný pilot pro manuální splice PS1.

## Co patří do registru

Registr přijímá pouze missense P/LP reference v kanonickém ENIGMA transkriptu,
které mají dohledatelný klasifikační původ. Přípustné klasifikační základy jsou:

1. P/LP reference z oficiální ENIGMA Supplementary Table 7 v1.2;
2. verzovaná oficiální ENIGMA/ClinGen VCEP assertion, například z ClinGen
   Evidence Repository, pokud ještě není v ST7;
3. lokální úplná reklasifikace podle uvedené verze ENIGMA VCEP pravidel.

Lokální reklasifikace musí být označena
`locally_recurated_under_enigma_vcep`. Nesmí se vydávat za oficiální expert
panel assertion. Musí obsahovat identifikátor posuzovatele, datum, verzi
pravidel a identifikátor úplného evidenčního záznamu.

Nové oficiální nebo lokálně reklasifikované reference mimo ST7 se zapisují do
`backend/data/ps1_protein_reference_extensions.json`. Generátor je sloučí se
ST7 a vytvoří jediný runtime registr. Prázdný extension soubor znamená, že
aktuální registr obsahuje pouze oficiální ST7 reference.

Samotný záznam v ClinVar bez ENIGMA/ClinGen expert-panel assertion, CANVarUK,
BRCA Exchange, jednotlivá publikace nebo výpočetní predikce nestačí k vytvoření
`eligible` reference. Tyto zdroje lze použít k nalezení kandidáta nebo jako
podklad úplné reklasifikace.

## Odkud se berou jednotlivá pole

| Informace v registru | Zdroj |
|---|---|
| gen, c. a p. notace, P/LP třída a původní klasifikační zdroj | ST7; u rozšíření oficiální VCEP assertion nebo lokální evidenční záznam |
| referenční transkript | BRCA1 `NM_007294.4`, BRCA2 `NM_000059.4` podle ENIGMA v1.2 |
| normalizovaná missense substituce | kanonická p. notace ověřená proti referenčnímu transkriptu |
| podklad proteinového mechanismu | PS3 funkční evidence z Table 9, nebo u patogenní missense reference doložená absence predikovaného a potvrzeného splice efektu |
| SpliceAI reference | pro současné ST7 sestavení pole `spliceai_prediction` z oficiální Table 9 v1.2 |
| známá RNA/splice evidence | úplná ENIGMA Table 9 v1.2 a úplná Supplementary Table 2 v1.2 |
| stav `eligible`, `excluded`, `review_required` | deterministické rozhodnutí generátoru z typu varianty, SpliceAI a známé splice evidence |
| případná PS1 závislost klasifikace reference | oficiální assertion nebo úplný lokální evidenční záznam |
| provenance a checksumy | generátor registru ze všech použitých verzovaných vstupních souborů |

Z 35 současných `eligible` referencí má 30 v Table 9 PS3 Strong funkční
evidenci. U zbývajících pěti je proteinový mechanismus zaznamenán jako
patogenní missense reference bez predikovaného nebo potvrzeného splice efektu.
Přímý proteinový funkční test tedy není povinný pro každý záznam, ale použitý
mechanistický podklad musí být v registru explicitní.

Pro budoucí SpliceAI zdroj musí být uložena verze modelu, referenční genom,
transkript, vstupní varianta, datum a checksum vstupních dat. Živá API odpověď
bez takové provenance sama nestačí pro trvalý záznam `eligible`.

Pro hodnocenou variantu uvedenou v Table 9 se používá její zmrazená
VCEP-revidovaná hodnota `spliceai_prediction`. Má přednost před obecnou
SpliceAI cache nebo API. Obecná cache se použije pro varianty, které v Table 9
revidované nejsou.

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

Po úmyslné kurátorské změně schváleného záznamu lze jeho checksum přepočítat
výhradně explicitním příkazem:

```bash
python scripts/validate_ps1_protein_registry.py --write-checksums
```

Přepočítání checksumu samo nenahrazuje odborné schválení obsahu.

## Kruhová závislost

Pokud je známo, že klasifikace reference použila PS1, musí být uložena použitá
reference. Bez ní záznam nelze schválit. Validátor odmítá přímou závislost na
sobě a známé cykly. Není však nutné u každé klasifikace dokazovat, že PS1
nepoužila, pokud použití PS1 není známé.
