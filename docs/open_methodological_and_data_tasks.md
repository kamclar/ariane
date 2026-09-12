# Otevřené metodické a datové úkoly

Tento dokument obsahuje pouze body, které ještě vyžadují externí metodické
potvrzení, nový validovaný dataset nebo dokončení validační práce. Opravené
auditní nálezy sem nepatří.

## 1. PM5 PTC uvnitř posledního exonu

Table 4 rozlišuje v posledním exonu dvě větve podle polohy PTC, zatímco
Appendix D určuje sílu podle exonu nukleotidové změny a uvádí příklad PTC v
následujícím exonu. Není potvrzeno, zda se v hraniční větvi mají PVS1 a PM5 PTC
posuzovat odděleně. Dotaz byl odeslán Janě. Do obdržení odpovědi se současná
logika nesmí měnit odhadem.

## 2. Úplný seznam patogenních founder variant

ENIGMA neposkytuje úplný strojově čitelný katalog founder variant, pro které se
nesmějí použít BA1 a BS1. Současný seznam obsahuje pouze doložené záznamy.
Pozitivní nález BA1 a BS1 vyloučí. Nenalezení má stav `not_listed` a znamená
pouze negativní výsledek tohoto konkrétního screeningu. Nejde o obecné tvrzení
o founder statusu.
Rozšíření vyžaduje kanonickou HGVS notaci, referenční transkript, populaci,
tvrzení o patogenitě, stabilní zdroj, datum přístupu a checksum.

## 3. Nezávislý validační soubor

Externě se připravuje soubor očekávaných ENIGMA klasifikací. Po jeho získání je
třeba připnout verzi a původ, oddělit automatizovatelnou a manuální evidenci a
porovnat jednotlivá kritéria, síly, body, výslednou třídu a důvody rozdílů.

## 4. Data pro automatické proteinové PS1, vyřešeno 2026-09-07

Odborná metodická konzultace potvrdila, že P/LP zařazení v ENIGMA ST7 v1.2 lze
považovat za splnění klasifikačního požadavku reference pro proteinové PS1.
Registr proto obsahuje 40 záznamů `eligible` a 20 záznamů `excluded` kvůli
známému splice efektu. Způsobilá reference přidá body pouze při stejné
normalizované missense substituci z jiné nukleotidové změny, SpliceAI nejvýše
0,1 u obou variant a splnění všech zaznamenaných proteinových a splice
podmínek. Původní otázka a rozhodnutí jsou popsány v
[`ps1_reference_validation_request.md`](ps1_reference_validation_request.md).

## 5. Nekvantifikovaná RNA evidence, částečně vyřešeno 2026-09-12

Nekvantifikované výsledky ST2 nyní správně vedou do předvyplněné manuální
revize a samy nepřidávají PVS1 RNA. Automatická cesta byla doplněna pouze pro
přesné publikované assertions ClinGen ERepo, které ENIGMA BRCA1/2 VCEP vytvořil
podle specifikace v1.2 a ve kterých je konečná RNA síla výslovně určena.

Runtime používá checksumovaný lokální registr, ne živý fallback. První registr
je záměrně označen jako neúplný. Zbývá zavést pravidelný export nových ERepo
záznamů, odbornou kontrolu kandidátů a vydání nové verze registru s checksumem.
Do té doby nepřítomnost varianty v registru nesmí být interpretována jako
nepřítomnost RNA evidence.

## 6. PP4/BP5: rozdíl mezi combined LR a štítkem zdrojového tracku

U 40 variant se síla vypočtená z nezkráceného `combined LR` podle publikovaných
prahových hodnot ENIGMA VCEP v1.2 liší od štítku `ACMGcode` v UCSC ENIGMA
tracku. Jde o 1 variantu u hranice BP5 Very Strong `0,00285`, 27 variant u
hranice BP5 Strong `0,05`, 8 variant u hranice BP5 Moderate `0,23` a 4 varianty
u hranice PP4 Moderate `4,3`. U posledních čtyř je výpočet ARIANE o jeden
stupeň silnější než štítek zdroje. U ostatních 36 je výpočet ARIANE o jeden
stupeň slabší.

Rozdíly odpovídají dvěma zveřejněným sadám prahů. Specifikace ENIGMA VCEP
v1.2, vydaná 9. ledna 2025, uvádí pro PP4 a BP5 hranice `4,3`, `0,23`, `0,05`
a `0,00285`.
Zanti et al., publikováno 25. května 2025, používá obecné kalibrované hranice
`4,33`, `0,231`, `0,053` a `0,0029`. UCSC v srpnu 2026 nahradilo původní
case-control složku hodnotami Zanti, přepočítalo combined LR a současně uvádí,
že multifaktoriální likelihood track je sestaven z publikovaných studií a je
nezávislý na verzi specifikace.

Srpen 2026 tedy není datem vydání nové VCEP specifikace. Jde o datum
aktualizace UCSC datového tracku. Aktuální registr ClinGen byl zkontrolován
1. září 2026 a jako nejnovější schválenou specifikaci stále uvádí verzi 1.2,
schválenou 9. ledna 2025. Záznam stejné verze na Zenodu vznikl technicky
18. července 2026, ale nemění datum schválení ani obsah pravidel.

Příkladem je BRCA2 `c.7805+34T>G`. Track uvádí combined LR `0,00286` a štítek
BP5 Very Strong. Hodnota splňuje hranici Zanti `<= 0,0029`, ale nesplňuje
hranici VCEP v1.2 `<= 0,00285`. Podle VCEP v1.2 proto ARIANE přiřadí
BP5 Strong. Nejde o rozdílnou hodnotu LR, ale o rozdílnou hranici pro převod LR
na sílu kritéria.

Zanti je novější než VCEP v1.2, ale novější publikace sama nemění platnou
specifikaci VCEP. ARIANE proto používá numerické LR z aktuálního tracku a sílu
odvozuje z prahů připnuté politiky VCEP v1.2. Původní `ACMGcode`, použitá sada
prahů a vypočtený výsledek VCEP se zachovávají pro audit. Prahy ani body se
nemění podle štítku zdrojového datasetu.

Zbývající dotaz k metodickému potvrzení: Má se aktualizovaný combined LR z
UCSC tracku vždy převést na sílu pomocí doslovných prahů CSpec v1.2, i když
štítek tracku odpovídá později publikovaným hranicím Zanti? Do vydání nové
verze VCEP nebo získání výslovného metodického potvrzení ARIANE používá prahy
VCEP v1.2.

Zdroje: [ClinGen ENIGMA BRCA1/2 VCEP v1.2](https://cspec.genome.network/cspec/ui/svi/doc/GN092?version=1.2.0),
[Zanti et al. 2025](https://www.nature.com/articles/s41467-025-59979-6) a
[UCSC ENIGMA track](https://hgdownload.soe.ucsc.edu/hubs/enigma/enigma.html).

## 7. Lokální HGVS mapování pro další geny

Současný referenční balík podporuje lokální převod `c.` na `p.`, ale není
obecným zdrojem genomových souřadnic. Obsahuje pouze transkriptový alignment pro
GRCh38 a neobsahuje genomové sekvence. Produkční resolver proto nyní používá
checksumované lokální souřadnicové mapy registrované v manifestu.

Před rozšířením na větší počet genů je vhodné připravit lokální HGVS
provider s připnutými genomovými sekvencemi a transkriptovými alignmenty pro
GRCh37 i GRCh38. Nový provider musí jednoznačně kontrolovat referenční alelu,
normalizaci indelů, transkript a assembly. Nejednoznačný nebo neúplný výsledek
musí zůstat nedostupný. Po validaci lze provider připojit přes existující
souřadnicový manifest bez změny klasifikačního DAGu.

## Již opravené auditní body

- ST7 P/LP je přijatý klasifikační základ proteinové PS1 reference; body se
  přidělí až po splnění všech nezávislých identity, mechanismu a splice kontrol.
- Nekvantifikovaná ST2 sama nepřiděluje PVS1 RNA.
- Rozdílné vícenásobné BayesDel hodnoty se neslučují maximem.
- Malé indely do 50 bp jsou odděleny od strukturální větve Appendix G.
- Nenalezení v neúplném founder registru se nepovažuje za negativní výsledek.
- Coverage s nepotvrzenou kompatibilitou zůstává auditní údaj a nepřiděluje body.
- Povolené typy PP3, BP4, BP7 a BP1 odpovídají Figure 1A. Nepotvrzený in-frame
  indel ani intronická pozice bez ověřitelného vyloučení `+/-1,2` nevstupují do
  bodování. Hranice `+7/-21` a interakce s RNA a funkční evidencí mají úplnou
  regresní matici.

## 8. Současné použití PS3 nebo BS3 a PVS1 RNA, čeká se na VCEP

U variant na posledním nukleotidu exonu, například `BRCA1 c.5074G>A/C` nebo
`BRCA2 c.7976G>A/C`, může Specifications Table 9 uvádět PS3 Strong z assay
hodnotícího mRNA i protein a Supplementary Table 2 současně obsahovat samostatný
RNA výsledek. Není potvrzeno, kdy lze PS3 nebo BS3 a PVS1 RNA použít současně,
ani jak ENIGMA požaduje doložit nezávislost studií a assayů, aby se stejný
mechanismus nezapočítal dvakrát.

Dotaz pro odbornou konzultaci:

> U variant na poslednim nukleotidu exonu, napriklad BRCA1 c.5074G>A/C nebo
> BRCA2 c.7976G>A/C, uvadi Specifications Table 9 PS3 Strong z assay
> zachycujiciho mRNA i protein a ST2 samostatnou RNA evidenci. Lze PS3 Strong a
> PVS1 RNA zapocitat soucasne, pokud pochazeji z odlisnych studii nebo assayu?
> Pokud kombinovany funkcni assay zachycuje oba mechanismy, ma se PVS1 RNA
> pouzit pouze ze samostatne nezavisle mRNA evidence? Jakou minimalni
> dokumentaci nezavislosti ENIGMA vyzaduje, aby nedoslo k dvojimu zapocteni?

### Odpověď Jany, 9. září 2026

Jana uvedla:

- u `BRCA1 c.5074G>A/C` je PVS1 určeno pouze pro mRNA evidenci, například RNA
  z pacientského materiálu nebo minigene assay;
- PS3 je u těchto variant použitelné podle Findlay et al., protože assay
  prokázal ztrátu buněčné funkce BRCA1;
- pokud by assay prokázal pouze aberantní sestřih, použila by pouze PVS1;
- u `BRCA2 c.7976G>A/C` platí stejný postup. PS3 vychází z Biswas et al., kde
  byl hodnocen rescue lethality u BRCA2-null myší;
- PVS1 u `BRCA1 c.5074G>A` vychází z podkladů ST2 a ST3;
- současné použití PVS1 a PS3 považuje za biologicky neintuitivní, ale za
  odpovídající současným podkladům, pokud PS3 dokládá funkční dopad produktu
  aberantního sestřihu. Otázku předá VCEP k potvrzení.

Z odpovědi pro ARIANE vyplývá následující pracovní postup:

- PVS1 RNA se používá pouze pro mRNA evidenci, například RNA z pacientského
  materiálu nebo minigene assay;
- pokud assay pouze prokáže aberrantní sestřih, použije se PVS1 RNA bez PS3;
- PS3 lze zachovat, pokud další assay prokáže funkční dopad nad rámec pouhé
  detekce sestřihu, například ztrátu buněčné funkce nebo rescue lethality;
- u `BRCA1 c.5074G>A/C` tento funkční podklad představuje Findlay et al.;
- u `BRCA2 c.7976G>A/C` jej představuje Biswas et al.

Jana současně uvedla, že biologická interpretace kombinace není intuitivní a
zařadí otázku k vyjasnění směrem k VCEP. Bod proto zůstává otevřený. Do odpovědi
VCEP se přesný ST2 záznam pouze předvyplní k odborné revizi a automaticky
nepřidělí PVS1 RNA ani body. Explicitní PS3 z Specifications Table 9 zůstává
zachováno. Pokud odborník přijme PVS1 RNA současně s PS3, ARIANE ponechá obě
kritéria a zobrazí povinné upozornění na rozsah assay, funkční výsledek a
nezávislost evidence.

Do předvyplnění se přenášejí také všechny odpovídající publikace a assay údaje
ze Supplementary Table 3. ST3 slouží k dohledání podkladů a sama nepřiděluje
PVS1 RNA ani jeho sílu.

## 9. Proteinový delins s terminačním kodonem

V připnutém indelovém snapshotu je 22 coding DNA `delins` variant, jejichž
normalizovaný proteinový následek je proteinový `delins` obsahující `Ter`,
například `BRCA1 c.3789_3790delinsTT`
`p.(Leu1263_Lys1264delinsPheTer)`. HGVS takový zápis připouští, pokud je
terminační kodon součástí vložené proteinové sekvence. Nejde o jednoduchou
nonsense substituci.

ENIGMA Table 4 popisuje PTC váhy pro protein termination variants, ale její
přehled automatických větví uvádí nonsense a frameshift. Není výslovně určeno,
zda se proteinový delins s `Ter` má posoudit v PVS1 a PM5 PTC větvi, nebo jako
in-frame delins ve Figure 1A. ARIANE tyto varianty bez metodického potvrzení
nepřeklápí do PTC větve. Bod vyžaduje konzultaci s VCEP.

## 10. Strojově ověřitelná verze služby SpliceAI

ARIANE připíná scoring profil, kontroluje GRCh38, `distance=10000`, `mask=0`,
úplnost všech skórovacích komponent a přesnou verzi referenčního transkriptu.
Veřejná odpověď Broad SpliceAI Lookup však neobsahuje neměnný identifikátor
verze modelu, anotace ani image digestu. Z odpovědi proto nelze strojově
prokázat, že provozovatel veřejného endpointu nezměnil interní release při
zachování stejné URL.

Technická část byla vyřešena 2026-09-11. ARIANE používá vlastní on-demand službu
z image připnutého úplným SHA-256 digestem. Varianta se skóruje až při požadavku.
Instalace nejprve ověří přesný referenční transkript, parametry Appendix J a
všechna delta, REF a ALT skóre proti verzovanému validačnímu případu. Teprve pak
přepne ARIANE na lokální endpoint. Změna image vyžaduje novou verzi profilu a
oddělený prostor runtime cache. Veřejný Broad endpoint není záložní zdroj.

Zůstává licenční otázka před případným zpřístupněním komerčním klinickým
laboratořím. Wrapper SpliceAI Lookup má licenci MIT. Použitý commit SpliceAI má
zdrojový kód pod GPLv3 a modelové váhy pod CC BY-NC 4.0. Současný bezplatný
akademický vývoj a interní hodnocení odpovídají nekomerčnímu účelu. Před použitím
v placené diagnostické službě je potřeba písemné vyjasnění s Illumina.

Podrobnosti implementace jsou v
[`implementation_and_data_sources.md`](implementation_and_data_sources.md).
