# Kontrola kritérií uvedených v tutorialových přepisech

Datum kontroly: 2026-08-04

Tento audit vychází přímo z textových přepisů uložených v
`F:/UOCHB/Enigma/Educational`. Sloupec „tutorial“ zachycuje to, co řečníci
navrhují při použití ENIGMA/VCEP pravidel. Nejde o přepis klasifikace z XLS.

| Run | Varianta | Kritéria navržená tutorialem | Třída tutorialu | Čas v přepisu |
|---|---|---|---:|---|
| 14 | BRCA1 c.509G>A p.(Arg170Gln) | BS1 Supporting; BP1 Strong; BS3 Strong | 1 | 24:08 až 30:49 |
| 14 | BRCA1 c.1534C>T p.(Leu512Phe) | BS1 Supporting; BP1 Strong; BS3 Strong; BP5 Strong | 1 | 31:18 až 41:07 |
| 14 | BRCA1 c.3668_3671dup p.(Cys1225fs) | PVS1 Very Strong; PM2 Supporting; PM5 PTC Strong | 5 | 41:45 až 48:00 |
| 14 | BRCA2 c.9097del p.(Thr3033fs) | PVS1 Very Strong; PM2 Supporting; PM5 PTC Strong | 5 | 48:56 až 50:12 |
| 15 | BRCA1 c.5551_5552insT p.(Asp1851ValfsTer29) | PVS1 Very Strong; PM5 PTC nepoužít; PM2 nepoužít podle ENIGMA | 3 | 11:20 až 22:00 |
| 15 | BRCA2 delece exonu 10 | PM2 Supporting; BS3 Moderate; PVS1 nepoužít | 3 | 22:17 až 29:39 |
| 15 | BRCA2 c.6147_6149del p.(Val2050del) | BP1 Strong | 2 | 29:49 až 37:00 |
| 13 | BRCA1 c.3891_3893del p.(Ser1298del) | BS3 Strong; BP1 Strong; BP5 Supporting | 1 | 17:36 až 28:22 |
| 13 | BRCA1 c.4185G>A p.(Gln1395=) | PVS1 RNA Strong; PM2 Supporting; PP4 Very Strong | 5 | 28:22 až 44:14 |
| 13 | BRCA1 c.628C>T p.(Gln210Ter) | PM2 Supporting; PVS1 nepoužít | 3 | 44:14 až 49:40 |
| 13 | BRCA2 c.8953+2T>C p.(?) | PM2 Supporting; PVS1 nepoužít; PP3 nepoužít podle ENIGMA | 3 | 49:40 až 57:15 |

## Současný výstup ARIANE

Externí ClinVar a ClinGen odpovědi byly při kontrolním běhu vypnuté. Klasifikace
proto vychází z lokálních verzovaných dat a stejných automatizovaných pravidel
jako produkční klasifikační cesta.

| Varianta | ARIANE kritéria | ARIANE třída | Porovnání s tutorialem |
|---|---|---:|---|
| BRCA1 c.509G>A | BS1 Supporting; BS3 Strong; BP1 Strong; clinical LR review required | 1 | Třída sedí. Multifaktoriální a Zanti case-control LR se automaticky nenásobí, dokud není doložena nezávislost zdrojových balíků. |
| BRCA1 c.1534C>T | BS1 Strong; BS3 Strong; BP1 Strong; clinical LR review required | 1 | Třída sedí. Kandidátní BP5 se automaticky neboduje kvůli nedoložené nezávislosti zdrojových balíků. |
| BRCA1 c.3668_3671dup | PVS1 Very Strong; PM5 PTC Strong | 5 | Třída sedí; tutorial používá také PM2 Supporting. ARIANE PM2 pro malé indely nepoužívá. |
| BRCA2 c.9097del | PVS1 Very Strong; PM5 PTC Strong | 5 | Třída sedí; tutorial používá také PM2 Supporting. |
| BRCA1 c.5551_5552insT | PVS1 Very Strong; PM5 PTC Strong | 5 | Neshoda s historickým tutorialem. Současná implementace ENIGMA v1.2 vybírá PVS1 a PM5 ze stejného řádku Table 4 podle exonu nukleotidové změny. |
| BRCA2 delece exonu 10 | PM2 Supporting; PVS1 N/A zobrazeno jako vyloučené | 3 | Třída se shoduje, kritéria ne. PM2 vzniká obecným Appendix G grafem nad Table 4 exony a úplným gnomAD-SV, nikoli variantovým záznamem. Tutorial používá BS3 Moderate, ale delece není v ENIGMA Table 9, proto ARIANE BS3 nepřidělí. |
| BRCA2 c.6147_6149del | BS1 Supporting; BP1 Strong | 2 | Třída sedí; BS1 Supporting je v ARIANE navíc. |
| BRCA1 c.3891_3893del | BS3 Strong; BP1 Strong; clinical LR review required | 1 | Třída sedí. Kandidátní klinický LR se neboduje kvůli nedoložené nezávislosti zdrojových balíků. |
| BRCA1 c.4185G>A | PVS1 RNA Strong; PM2 Supporting; PP4 Strong | 4 | PVS1 RNA se obecně odvodilo z přesného ST2 řádku, delece exonu 12, Table 4 a kvalitativní větve Appendix E. PP3 bylo potlačeno jako slabší evidence stejného splice mechanismu. Proti tutorialu zůstává rozdíl PP4 Strong versus Very Strong. |
| BRCA1 c.628C>T | PM2 Supporting | 3 | Shoda. Table 4 uvádí PVS1 N/A. |
| BRCA2 c.8953+2T>C | PM2 Supporting | 3 | Shoda. Table 4 uvádí PVS1 N/A a ENIGMA tutorial PP3 nepoužívá. |

## Důležité poznámky k přepisu

1. Automatický přepis obsahuje chyby v názvech variant a zkratkách. Například
   u `p.Arg170Gln` řečník podle přepisu říká změnu na glycin. Identita varianty
   byla převzata ze zadání tutorialu, kritéria z navazujícího výkladu.
2. U `c.4185G>A` tutorial v čase 30:33 výslovně odkazuje na Parsons et al. 2019
   a automatický přepis zachytil LR jen jako „8.6“. Parsonsovy komponenty
   segregation `39,31669` a pathology `13,9129` však dávají LR `547,009`, tedy
   PP4 Very Strong. Přepis proto zřejmě zkomolil nebo vynechal začátek vyslovené
   hodnoty. Současný oficiální ENIGMA BRCAmfa track přidává nezávislou komponentu
   personal/family history `0,59996` z Li et al. 2020. Výsledný combined LR je
   `328,184`, tedy PP4 Strong. To vysvětluje rozdíl mezi tutorialem a současným
   výstupem ARIANE i HECTORu.
3. U `c.3891_3893del` tutorial výslovně používá posterior probability 0,368 pro
   BP5 Supporting a odmítá tehdy dostupný LR 28 pro PP4 Strong. ARIANE eviduje
   komponenty Parsons 2019, Caputo 2021 a Zanti 2025. Jejich kandidátní součin
   `0,0289608` se automaticky neboduje, protože nezávislost zdrojových balíků
   není doložena ve verzované matici překryvu.
4. U `c.5551_5552insT` tutorial výslovně vysvětluje, že PM5 PTC nelze použít podle polohy předpokládaného terminačního kodonu. Appendix D ENIGMA v1.2 ale uvádí, že sílu PM5 určuje exon nukleotidové změny, a u frameshiftu se stop kodonem v následujícím exonu používá pravidlo původního exonu. ARIANE proto nyní páruje PVS1 a PM5 ze stejného řádku Table 4. Tento rozdíl je v auditu zachován jako rozdíl historického tutorialu a současné implementace v1.2.
   První změněná aminokyselina je 1851 a `fsTer29` umísťuje předpokládaný stop
   za hranici 1854. ARIANE obě polohy uvádí v auditu, ale pro výběr páru PVS1 a
   PM5 používá řádek Table 4 pro exon nukleotidové změny.
5. `BRCA1 c.3247A>C p.(Met1083Leu)` a `BRCA1 c.5217T>A p.(Asp1739Glu)` nejsou
   v těchto tutorialových přepisech vysvětleny. Jejich kritéria nelze vydávat za
   tutorialový výsledek; pocházejí pouze z XLS nebo z jiné srovnávací sady.
