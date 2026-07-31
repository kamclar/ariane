# Matice interakcí funkční, RNA a predikční evidence

## Rozsah

Tato matice řídí mechanism-aware deduplikaci v ARIANE. Není založena na obecném příznaku přítomnosti funkční evidence. Rozlišuje proteinovou funkci, experimentálně hodnocený dopad na mRNA, predikovaný splice efekt, predikovaný proteinový efekt a doménový kontext.

Primární zdroje:

- [ENIGMA BRCA1/2 Specifications v1.2](https://cspec.genome.network/cspec/File/id/11e62fec-23b0-4a3e-b2df-751855301746/data), Figure 1A, Figure 1B a Figure 1C
- [ENIGMA BRCA1/2 Appendix v1.2](https://cspec.genome.network/cspec/File/id/9e6119dc-90b9-42b5-a3b7-1a2eb28b1b12/data), Appendix E
- ENIGMA Specifications Table 9 pro kalibrovaná PS3 a BS3

## Mechanismy kritérií

| Kritérium | Mechanismus nebo komponenty |
| --- | --- |
| PP3 ze SpliceAI | predikovaný splice efekt |
| PP3 z BayesDel | predikovaný proteinový efekt v doméně |
| BP4 pro missense/in-frame | negativní proteinová i splice predikce v doméně |
| BP4 pro silent/intronic | negativní splice predikce |
| BP7 Supporting | negativní splice predikce plus nízký prior silent nebo hluboké intronické varianty |
| BP1 Strong | poloha mimo funkční doménu plus negativní splice predikce |
| PS3/BS3 | kalibrovaný dopad na proteinovou funkci, s povinným zohledněním splice kontextu |
| PVS1 (RNA) | experimentálně potvrzený škodlivý dopad na mRNA |
| BP7 Strong (RNA) | experimentálně potvrzená absence škodlivého dopadu na mRNA |

BP1 není čistě proteinové kritérium. Obsahuje doménovou komponentu i podmínku SpliceAI nejvýše 0,1.

## Rozhodovací matice

| Silnější nebo přidaná evidence | Slabší nebo souběžná evidence | Akce | Zdroj |
| --- | --- | --- | --- |
| automatické PVS1 | PP3 | PP3 nezapočítat | Figure 1A a 1B |
| přijaté PVS1 (RNA) | PP3, BP4, BP7 nebo BP1 | bioinformatické kódy nezapočítat | Figure 1B: replace bioinformatic codes |
| přijaté PVS1 (RNA) | proteinové PS1 nebo predikční PS1 splicing | PS1 nezapočítat, pokud popisuje nahrazený proteinový nebo splice předpoklad | Figure 1B |
| přijaté PVS1 (RNA) | PS3 nebo BS3 | zachovat, ale vyžádat revizi rozsahu a nezávislosti assay | Figure 1C a Appendix E |
| přijaté BP7 Strong (RNA) | BP7 Supporting | BP7 Supporting nahradit silnějším BP7 Strong (RNA) | Figure 1B a Appendix E |
| přijaté BP7 Strong (RNA) | BP1 nebo BP4 | zachovat, pokud jsou jinak použitelné | Figure 1B |
| přijaté BP7 Strong (RNA) | splice PP3 | zachovat obě položky, označit konflikt a vyžádat revizi | Figure 1B |
| přijaté BP7 Strong (RNA) | PS3 nebo BS3 | zachovat a vyžádat revizi podmínek Figure 1B/1C | Figure 1B, Figure 1C a Appendix E |
| PS3 nebo BS3 bez PVS1 | relevantní PP3, BP4, BP7 nebo BP1 | zachovat relevantní bioinformatické kódy | Figure 1C |
| PS3/BS3 a bioinformatický kód v opačném směru | všechny dotčené kódy | zachovat, označit konflikt a vyžádat expertní revizi | Figure 1C |
| nejasný mRNA výsledek | bioinformatické kódy | neměnit, uvést nejasný RNA výsledek v popisu | Figure 1B |

## Fail-closed pravidla

- Mechanismus se nesmí odhadnout pouze z existence PS3 nebo BS3.
- Chybějící informace o rozsahu assay nesmí automaticky odstranit kritérium.
- Deduplikace se provede pouze tam, kde Figure 1B výslovně stanoví nahrazení nebo upgrade.
- Při konfliktu experimentální a predikční evidence se položky zachovají, pokud ENIGMA flowchart nestanoví nahrazení.
- Každé odstraněné kritérium musí být uvedeno ve strukturovaném `evidence_interactions` včetně zachovaného kritéria, mechanismu, důvodu a odkazu.

## Uživatelský výstup

Veřejný i amended výsledek může obsahovat rozbalovací část `Evidence interaction warnings`.

Stavy:

| Stav | Význam |
| --- | --- |
| `info` | kódy jsou podle flowchartu zachovány a uživateli se vysvětluje jejich vztah |
| `deduplicated` | slabší důkaz nebyl započítán, protože jej nahradil silnější důkaz stejného mechanismu nebo následku |
| `review_required` | kombinace může být přípustná, ale vyžaduje kontrolu assay a mechanismu |
| `conflict` | důkazy pro stejný nebo související mechanismus ukazují opačným směrem |

