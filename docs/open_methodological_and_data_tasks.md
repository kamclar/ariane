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

## 2. gnomAD coverage pro PM2

ENIGMA vyžaduje dostatečnou mean depth v oblasti kolem varianty, ale neurčuje
šířku tohoto okolí. ARIANE nyní používá přesně genomový rozsah alely `REF` a
uvádí jej jako `coverage_scope: variant_reference_span`. Tato hodnota zůstává
dostupná pro audit, ale není považována za potvrzenou metodu ENIGMA a sama
nemůže automaticky založit PM2.

Současně jsou frekvence z gnomAD v3.1.2 kombinovány s veřejně dostupnou coverage
z release 3.0.1. Runtime tuto dvojici označuje jako
`unresolved_release_mismatch`. Hloubku zobrazí, ale nepoužije ji k automatickému
BA1, BS1 ani PM2. Je třeba získat metodické potvrzení přijatelné šířky okolí a
ověřit, zda je tato kombinace verzí přípustná. Případná změna musí mít novou
verzi datové politiky a regresní testy.

## 3. Úplný seznam patogenních founder variant

ENIGMA neposkytuje úplný strojově čitelný katalog founder variant, pro které se
nesmějí použít BA1 a BS1. Současný seznam obsahuje pouze doložené záznamy.
Nenalezení varianty proto vrací stav `unresolved`, nikoliv negativní závěr.
BA1 nebo BS1 lze použít jen po explicitním výsledku `reviewed_not_found`.
Rozšíření vyžaduje kanonickou HGVS notaci, referenční transkript, populaci,
tvrzení o patogenitě, stabilní zdroj, datum přístupu a checksum.

## 4. Nezávislý validační soubor

Externě se připravuje soubor očekávaných ENIGMA klasifikací. Po jeho získání je
třeba připnout verzi a původ, oddělit automatizovatelnou a manuální evidenci a
porovnat jednotlivá kritéria, síly, body, výslednou třídu a důvody rozdílů.

## 5. Data pro automatické proteinové PS1

ST7 je bezpečně používána jen jako zdroj kandidátů. Automatické PS1 bude možné
rozšířit, až bude pro konkrétní reference dostupná samostatně ověřená
ENIGMA/ClinGen VCEP assertion nebo úplná lokální reklasifikace podle uvedené
verze VCEP pravidel. Samotný záznam ST7 nestačí.

## 6. Nekvantifikovaná RNA evidence

Nekvantifikované výsledky ST2 nyní správně vedou do předvyplněné manuální
revize a samy nepřidávají PVS1 RNA. Další automatizace by vyžadovala
strukturovaný zdroj, který výslovně zachytí konsenzuální kurátorské zařazení do
větve Appendix E Table 9 a výslednou sílu. Bez takového zdroje zůstává tato část
manuální.

## 7. PP4/BP5: rozdíl mezi combined LR a štítkem zdrojového tracku

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

## 8. Lokální HGVS mapování pro další geny

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

- ST7 sama nevytváří automaticky způsobilou proteinovou PS1 referenci.
- Nekvantifikovaná ST2 sama nepřiděluje PVS1 RNA.
- Rozdílné vícenásobné BayesDel hodnoty se neslučují maximem.
- Malé indely do 50 bp jsou odděleny od strukturální větve Appendix G.
- Nenalezení v neúplném founder registru se nepovažuje za negativní výsledek.
- Coverage s nepotvrzenou kompatibilitou zůstává auditní údaj a nepřiděluje body.

Podrobnosti implementace jsou v
[`implementation_and_data_sources.md`](implementation_and_data_sources.md).
