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

## Již opravené auditní body

- ST7 sama nevytváří automaticky způsobilou proteinovou PS1 referenci.
- Nekvantifikovaná ST2 sama nepřiděluje PVS1 RNA.
- Rozdílné vícenásobné BayesDel hodnoty se neslučují maximem.
- Malé indely do 50 bp jsou odděleny od strukturální větve Appendix G.
- Nenalezení v neúplném founder registru se nepovažuje za negativní výsledek.
- Coverage s nepotvrzenou kompatibilitou zůstává auditní údaj a nepřiděluje body.

Podrobnosti implementace jsou v
[`implementation_and_data_sources.md`](implementation_and_data_sources.md).
