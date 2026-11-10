# Otevřené metodické a datové úkoly

Tento dokument obsahuje pouze body, které ještě vyžadují externí metodické
potvrzení, nový validovaný dataset nebo dokončení validační práce. Opravené
auditní nálezy sem nepatří.

## 1. PM5 PTC uvnitř posledního exonu, vyřešeno 2026-10-15

ENIGMA VCEP potvrdil, že v hraniční větvi posledního exonu musí první změněná
aminokyselina i nově zavedený terminační kodon ležet na kritické hranici nebo
před ní. PVS1 a PM5 PTC se posuzují společně podle jednoho řádku Table 4. Pokud
jedna z poloh podmínku nesplní nebo terminační pozici nelze určit, nepoužije se
ani jeden kód. Pokud bez terminační pozice nelze vybrat větev, ARIANE nevrátí
automatickou klasifikaci a vyžádá odbornou revizi. Implementace a regresní
případy jsou popsány v hlavním implementačním dokumentu.

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

## 5. RNA evidence, automatizovatelná část vyřešena 2026-11-10

Primární automatická cesta používá konečný alelový kód označený `(RNA)` přímo
z ENIGMA Specifications Table 4. Jde o výsledek VCEP po zohlednění pozorovaných
RNA produktů, nikoli o nový odhad ARIANE. Aktuální tabulka obsahuje 144 přesných
pravidel: 142 přiděluje PVS1 RNA v uvedené síle a dvě určují N/A. Přesné Table 4
rozhodnutí nahrazuje predikční splice větev; SpliceAI se proto nevyžaduje ani u
N/A záznamu, který nepřidává body.

U alel bez konečného Table 4 RNA kódu mohou přesné ST2 řádky v
ENIGMA-kurátorované kategorii pacientské mRNA bez alelově
specifické kvantifikace s aberantním transkriptem odpovídajícím ztrátě funkce
přidělit PVS1 RNA. Automatická větev navíc vyžaduje propojenou ST3
referenci a jednoznačné mapování delece celého exonu na Table 4. Síla se odvodí
z kvalitativní apparent near-complete větve Appendix E. Komplexní a částečné
transkriptové výsledky zůstávají v předvyplněné manuální revizi bez bodů pouze
tehdy, pokud pro danou alelu není konečný Table 4 nebo schválený ERepo výsledek.

Další automatická cesta používá přesné publikované assertions ClinGen ERepo,
které ENIGMA BRCA1/2 VCEP vytvořil podle specifikace v1.2 a ve kterých je
konečná RNA síla výslovně určena.

Synchronizace z 9. listopadu 2026 doplnila dvě takové aktuální v1.2 assertions:
`BRCA2 c.7976+5G>T` a `BRCA2 c.8486A>G`, obě s `PVS1 (RNA) Very Strong`.
Konzistenční test nyní vyžaduje, aby každá aktuální v1.2 ERepo assertion s
explicitní konečnou RNA vahou byla přítomna i ve schváleném runtime registru.
Záznamy bez výslovné finální síly nadále zůstávají pouze v předvyplněné odborné
revizi; jejich neautomatizování není implementační chyba.

Runtime používá checksumovaný lokální registr, ne živý fallback. První registr
je záměrně označen jako neúplný. Týdenní hlídač nyní porovnává ERepo se
zaznamenaným stavem zdroje. Při změně vytvoří nebo aktualizuje GitHub issue s
rozdílem. Samostatně spuštěný workflow umí stáhnout kandidátní snapshot,
porovnat záznamy a kritéria, spustit úplné testy a připravit draft pull request.
Nový registr se však neaktivuje automaticky: kandidáty musí odborně zkontrolovat
člověk a změna vstoupí do ARIANE až po schválení a sloučení pull requestu.
Nepřítomnost varianty v registru nadále nesmí být interpretována jako
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

## Stav automatizace aktualizace zdrojů

ARIANE odděluje zjištění změny, přípravu kandidátních dat, odborné schválení a
produkční nasazení. Novější vzdálený zdroj proto nikdy sám nepřepíše aktivní
snapshot ani klasifikační politiku.

- Týdenní source watch kontroluje ERepo, UCSC ENIGMA clinical LR, ClinGen CSpec
  a dostupnost odpovídajících produktů gnomAD. Změna zdroje je review událost,
  nikoli chyba CI. Vznikne deduplikované GitHub issue s identitou zdroje,
  zjištěným rozdílem a navrženým dalším krokem. Skutečné technické selhání
  kontroly zůstává červenou chybou workflow.
- Pro ERepo lze připravit kandidátní registr včetně rozdílu záznamů a kritérií.
  Pro UCSC clinical LR lze připravit kandidátní snapshot včetně kontroly
  checksumu nástroje `bigBedToBed` a rozdílu hodnot LR, PP4/BP5 a směru
  evidence. Obě cesty spustí úplné testy a vytvoří pouze draft pull request.
- Změna CSpec vytvoří metodický report s verzí, revizí, datem a stavem vydání.
  Neimplementuje automaticky nová pravidla, protože změnu musí nejprve odborně
  přeložit a validovat člověk.
- U gnomAD se sleduje identita aktivních Hail objektů a dostupnost stejného
  produktu v novější verzi. Automatický přechod ani přestavba snapshotu se
  neprovádí: je nutné posoudit release, kohortu, coverage, schéma a kompatibilitu
  populační politiky.
- SpliceAI profil a kurátorované lokální registry mají další ruční kontrolu
  naplánovanou na 1. února 2027. Referenční transkripty a souřadnicové zdroje na
  1. května 2027.

Aktivní data se změní teprve po odborném schválení, sloučení verzovaného pull
requestu a samostatném nasazení. Přesný provozní postup a odpovědnosti jsou v
[`source_freshness_and_update_policy.md`](source_freshness_and_update_policy.md);
technické vazby zdrojů jsou v
[`implementation_and_data_sources.md`](implementation_and_data_sources.md).

## Již opravené auditní body

- U smíšené patogenní a benigní evidence se bodový systém vyhodnotí, ale
  nepublikuje ne-VUS třídu, pokud ji podporuje jen jeden kód v daném směru.
  Například `PVS1 Very Strong + BP7 Supporting` zůstává VUS a vyžaduje odbornou
  diskusi. Druhý nezávislý patogenní kód dovolí použít vypočtenou bodovou třídu.
- BP1 a BP7 normalizují token typu varianty na hranici pravidla stejnou sdílenou
  funkcí jako PP3/BP4. Velikost písmen ani okolní mezery proto nemohou změnit
  použitelnost kritéria, sílu, body nebo větev auditního rozhodovacího stromu.
- Rozpoznaný typ varianty se již nepovažuje automaticky za klasifikovatelný.
  Přesné `5utr` a `3utr` SNV mají omezenou větev pro populační kritéria a
  přesnou variantově specifickou kurátorovanou evidenci; samotná UTR poloha kód
  nevytváří. UTR indely končí jako `unsupported_utr_edit` a `stop_lost` jako
  `unsupported_variant_type`. Neurčené `unknown`, `delins`, `deletion`,
  `insertion` a `duplication` končí jako `protein_consequence_unresolved`.
  Nepodporované větve nevracejí VUS a výsledek se neukládá do cache.
- ST7 P/LP je přijatý klasifikační základ proteinové PS1 reference; body se
  přidělí až po splnění všech nezávislých identity, mechanismu a splice kontrol.
- Nekvantifikovaná ST2 přiděluje PVS1 RNA jen v přesné standardizované LoF
  kategorii s propojeným ST3 zdrojem a jednoznačným Table 4 důsledkem.
- Rozdílné vícenásobné BayesDel hodnoty se neslučují maximem.
- Malé indely do 50 bp jsou odděleny od strukturální větve Appendix G.
- Nenalezení v neúplném founder registru se nepovažuje za negativní výsledek.
- Coverage s nepotvrzenou kompatibilitou zůstává auditní údaj a nepřiděluje body.
- Povolené typy PP3, BP4, BP7 a BP1 odpovídají Figure 1A. Nepotvrzený in-frame
  indel ani intronická pozice bez ověřitelného vyloučení `+/-1,2` nevstupují do
  bodování. Hranice `+7/-21` a interakce s RNA a funkční evidencí mají úplnou
  regresní matici.
- Coding DNA `delins` s ověřeným proteinovým následkem obsahujícím předčasný
  `Ter` se typuje jako `ptc_delins` a směruje do Table 4 větve PVS1/PM5 PTC,
  nikoli do in-frame Figure 1A. V připnutém indelovém snapshotu tuto cestu
  regresně pokrývá všech 22 takových záznamů.

## 8. Současné použití PS3 nebo BS3 a PVS1 RNA, vyřešeno

U variant na posledním nukleotidu exonu, například `BRCA1 c.5074G>A/C` nebo
`BRCA2 c.7976G>A/C`, může Specifications Table 9 uvádět PS3 Strong z assay
hodnotícího mRNA i protein a Supplementary Table 2 současně obsahovat samostatný
RNA výsledek. Původně nebylo potvrzeno, kdy lze PS3 nebo BS3 a PVS1 RNA použít
současně ani jak doložit nezávislost studií a assayů, aby se stejný mechanismus
nezapočítal dvakrát.

Tuto otázku následně uzavřelo metodické upřesnění ENIGMA VCEP z 15. října 2026.
Původní dotaz pro odbornou konzultaci:

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

ENIGMA VCEP v metodickém upřesnění z 15. října 2026 potvrdil, že oba kódy lze
zachovat, pokud dokládají odlišné příspěvky: PVS1 RNA pro aberantní transkript a
PS3 nebo BS3 pro samostatný funkční dopad nad rámec pouhého průkazu sestřihu.
Stejné experimentální pozorování se nesmí započítat dvakrát. ARIANE proto
explicitní PS3 nebo BS3 z Table 9 zachová a při současném přijetí PVS1 RNA
zobrazí povinné upozornění na rozsah assay, funkční výsledek, provenance a
nezávislost evidence. Toto rozhodnutí nemění konzervativní pravidlo, že samotné
předvyplnění formuláře nepřiděluje PVS1 RNA ani body.

Do předvyplnění se přenášejí také všechny odpovídající publikace a assay údaje
ze Supplementary Table 3. ST3 slouží k dohledání podkladů a sama nepřiděluje
PVS1 RNA ani jeho sílu.

Od verze 1.9.16 RNA review navíc porovnává strukturované Table 9 pole pro
funkční assay s polem publikovaných splice výsledků. Pokud se v obou rolích
objeví stejný PMID, zobrazí jej jako kontrolu provenance. Shoda PMID sama
kritérium nemaže, protože publikace může obsahovat odlišné experimenty. Před
současným přijetím PS3/BS3 a PVS1 RNA je však nutné přiřadit konkrétní assay a
pozorování ke každému kódu a vyloučit opakované použití stejné evidence. U
`BRCA2 c.7976G>A` je takto označen PMID `33293522`; automaticky se započítává jen
PS3 a předvyplnění PVS1 RNA odkazuje na tři jiné ST3 publikace.

## 9. Strojově ověřitelná verze služby SpliceAI

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

## 10. Historická Class 5 a současná klasifikace podle VCEP v1.2

Varianta `BRCA1 NM_007294.4:c.4185G>A p.(Gln1395=)` má v ClinVaru ENIGMA
expert-panel assertion `SCV001161546` s klasifikací Pathogenic, hodnocenou
18. června 2019 podle ENIGMA BRCA1/2 Classification Criteria z roku 2017.
Komentář assertion uvádí IARC Class 5 na základě posteriorní pravděpodobnosti
`0,999943` z multifaktoriální likelihood analýzy. Nejde o zveřejněný
kritérium-po-kritériu rozpis podle současné specifikace VCEP v1.2.

ARIANE podle aktuálně připnutých podkladů v1.2 použije `PVS1 RNA Strong`,
`PP4 Strong` a `PM2 Supporting`. Combined LR `328,184` nedosahuje hranice `350`
pro `PP4 Very Strong`. Dvě silná a jedno podpůrné patogenní kritérium splňují
Table 3 kombinaci pro Likely Pathogenic, nikoli pro Pathogenic. Historická
výsledná multifaktoriální Class 5 se nepřidává jako další kritérium, protože by
se tím znovu započítaly její podkladové složky.

Současně zveřejnila ENIGMA BRCA1/2 VCEP v roce 2026 podle v1.2 kuraci varianty
`BRCA1 c.4185G>C`. Tato kurace výslovně používá `c.4185G>A` jako variantu
klasifikovanou ENIGMA VCEP jako Pathogenic a jako referenci pro `PS1 Splicing`.
ENIGMA tedy historický patogenní závěr nadále přijímá, ale pro přesnou variantu
`c.4185G>A` nebyl v kontrolovaném ERepo snapshotu nalezen nový samostatný v1.2
rozpis kritérií.

Do vyjasnění se historická klasifikace zobrazuje výrazně a odděleně od výpočtu:

> Historical ENIGMA expert-panel classification: Pathogenic; current ARIANE
> v1.2 automated result: Likely Pathogenic.

Externí assertion nemění body ani automatický výsledek. Rozdíl vyžaduje odbornou
revizi a následující dotaz pro ENIGMA:

> BRCA1 NM_007294.4:c.4185G>A p.(Gln1395=) has an ENIGMA expert-panel ClinVar
> assertion, SCV001161546, evaluated on 18 June 2019 as Pathogenic, IARC Class 5,
> based on a posterior probability of 0.999943 from multifactorial likelihood
> analysis under the 2017 ENIGMA criteria. Applying the evidence currently
> available to ARIANE under BRCA1/2 VCEP v1.2 gives PVS1 RNA Strong, PP4 Strong
> with a combined LR of 328.184, and PM2 Supporting. This combination meets
> Likely Pathogenic, but not Pathogenic, under v1.2 Table 3. A 2026 v1.2 ERepo
> assertion for BRCA1 c.4185G>C still cites c.4185G>A as a Pathogenic PS1
> Splicing reference. Should the historical expert-panel Class 5 for c.4185G>A
> remain authoritative under v1.2 without a new criterion-by-criterion
> assertion, or should c.4185G>A be re-curated under v1.2? If it remains
> Pathogenic, which current v1.2 criteria and evidence combination should be
> recorded?

Zdroje: [ClinVar RCV000112286](https://www.ncbi.nlm.nih.gov/clinvar/RCV000112286/)
a [ClinGen ERepo c.4185G>C](https://erepo.clinicalgenome.org/evrepo/ui/classification/9ef5ebb6-b717-4061-8655-b7c3bf79b9a5).

## 11. Dvě podezřelé alely ve zdrojové Table 4

Audit původního workbooku Table 4 našel dva případy, které vyžadují potvrzení
ENIGMA:

- BRCA2 obsahuje vzor `c.6937+1T>`. Referenční sekvence má na pozici `+1` bázi
  G a na pozici `+2` bázi T. Pravděpodobně má jít o `c.6937+2T>`.
- Pro BRCA2 `c.8487+2T>` jsou uvedeny alternativy C a G, ale chybí A. Je třeba
  potvrdit, zda jde o záměrnou alelovou výjimku nebo vynechanou alternativu.

ARIANE tyto případy neopravuje vlastní výjimkou. Zachovává zdrojový workbook a
pro chybějící přesnou alelu postupuje fail-closed. Případná oprava se provede
až v nové schválené verzi zdroje nebo na základě výslovného stanoviska ENIGMA.

Dvojí zápis `BRCA1 c.81-1G>C` a `BRCA2 c.8487+2T>C` není považován za chybu.
Workbook u nich obsahuje obecný predikční řádek a přesnější RNA řádek. Runtime
použije RNA rozhodnutí, ale oba původní řádky uchovává pro audit.
