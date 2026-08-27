# Implementace kritérií, klasifikace a datové zdroje ARIANE

## 1. Rozsah aplikace

ARIANE provádí první automatizovaný průchod pravidly ClinGen ENIGMA BRCA1/2 VCEP v1.2 pro geny BRCA1 a BRCA2.

Použité referenční transkripty:

| Gen | RefSeq | Ensembl pro předpočítaná data |
| --- | --- | --- |
| BRCA1 | `NM_007294.4` | `ENST00000357654.9` |
| BRCA2 | `NM_000059.4` | `ENST00000380152.8` |

### 1.1 Verzovaný registr genů a VCEP politik

Autoritativní runtime konfigurace je v
`backend/data/gene_policy_manifest.json`. Její metadata a SHA-256 jsou v
`backend/data/gene_policy_manifest.metadata.json`. Aplikace při startu odmítne
chybějící pole, neplatné pořadí prahů, neznámou politiku, neaktivní konfiguraci
nebo rozdílný checksum.

Manifest pro každý aktivní gen obsahuje:

- referenční transkript a protein,
- identifikaci a verzi VCEP specifikace,
- runtime policy ID,
- genově specifické prahy BayesDel_noAF a PVS1 NMD hranici,
- funkční domény,
- úplný seznam použitelných automatických a manuálních pravidel,
- požadované datové zdroje pro jednotlivé rodiny pravidel.

Sdílená VCEP politika obsahuje také verzovaný seznam původních ACMG/AMP použití,
která specifikace označuje jako `Do not use`. Každá položka má kód, krátký důvod
a část zdrojového dokumentu. Tento seznam se nesmí zaměňovat s kritériem, které
je pouze nesplněné nebo nepoužitelné pro konkrétní variantu. Seznam `Do not use`
je vlastnost VCEP politiky, nikoliv výsledek klasifikace jedné varianty.

Výsledek klasifikace má samostatná pole pro aplikovaná kritéria, kritéria
vyloučená výjimkou a kritéria explicitně označená jako `not_applicable` ve
variantově specifické rozhodovací cestě. Rozhraní zobrazuje poslední skupinu
v rozbalovacím řádku `Not applicable to this variant`. Pouhá absence kódu ve
výsledku se na `not_applicable` nepřevádí. Nesplněný práh, nedostupný zdroj,
nejednoznačná evidence a výjimkou vyloučené kritérium zůstávají odlišné stavy.
Prvními automaticky zveřejněnými případy jsou PVS1 označené jako N/A přímo v
ENIGMA Table 4 a PM2, jehož nepoužitelnost pro daný typ alely určuje aktivní
VCEP politika.

Sdílená část VCEP politiky obsahuje prahy SpliceAI, PP4, BP5, BA1, BS1, PM2 a
hranice bodové klasifikace pro mixed evidence. Python neobsahuje náhradní
BayesDel prahy pro neznámý gen. Gen bez aktivního záznamu je odmítnut. Kritérium,
které není uvedeno v `applicable_rules`, se nezapočítá a ve výsledku se objeví
jako vyloučené s důvodem.

Zdrojové manifesty gnomAD, SpliceAI a referenčního sekvenčního balíku zůstávají
oddělené, protože popisují konkrétní datasety. Startovní kontrola ale porovnává
jejich transkripty, VCEP identitu a rozhodovací prahy s hlavním registrem.
Rozpor zastaví aplikaci. Zdrojový manifest proto nemůže tiše změnit klasifikační
politiku.

Po schválené změně politiky nebo prahů se zvýší `manifest_version` a obnoví
kontrolní metadata příkazem:

```powershell
.\venv\Scripts\python.exe scripts\update_gene_policy_manifest_metadata.py --write
```

Přidání dalšího genu začíná novým záznamem v tomto manifestu. Je nutné dodat
také referenční sekvence, zdrojové datasety a regresní případy uvedené v
`required_rule_data`. Samotné přidání symbolu nebo gnomAD intervalu gen
neaktivuje. Pokud je v panelu více různých VCEP politik, požadavek bez genu se
odmítne jako nejednoznačný.

Parser genového prefixu, seznam genů ve formuláři a validační API nejsou omezené
regexem na BRCA1/2. Čtou seznam aktivních genů z registru. Každá politika má
navíc povinný `implementation_profile`. Produkční DAG přijme pouze profil, pro
který existuje explicitní implementace. Nový VCEP profil proto nemůže omylem
použít ENIGMA BRCA pravidla. Chybějící implementace ukončí požadavek chybou.

Genový záznam obsahuje také kontrolní variantu pro startup validaci HGVS,
genově specifické odkazy na VCEP specifikaci, PVS1 decision assets a popisy
funkčních domén. Tyto hodnoty již nejsou větvené podle názvu genu v Pythonu.
Názvy BRCA zůstávají pouze u zdrojů, které jsou skutečně BRCA specifické,
například ENIGMA Tables 4/9, ST7, founder varianty a exonové CNV. Takový dataset
se pro jiný gen nesmí použít bez odpovídajícího provideru a validace provenance.

Automatický výsledek není úplnou expertní klasifikací. Kritéria PS4, PM3, PP1, BS2 a BS4 vyžadují klinická, rodinná nebo literární data a automaticky se nepřidávají. Aplikace pro ně podporuje oddělenou strukturovanou manuální revizi. Uživatel zadává podklady, ale nemůže zvolit sílu kritéria. Backend ji vždy odvodí z prahů ENIGMA BRCA1/2 VCEP v1.2. Nenulové `override_strength` API výslovně odmítá a podprahové podklady nezískají kritérium ani body.

Manuální vyhodnocení kontroluje také povinné stipulace CSpec v1.2. PS4 vyžaduje shodu země a etnicity případů a kontrol. PM3 a BS2 vyžadují ověření, že koexistující P/LP varianta byla klasifikována podle VCEP specifications. PM3 navíc vyžaduje potvrzení, že hodnocená varianta nesplňuje benigní populační kritérium. PP1 Very Strong vyžaduje zaznamenaný predikovaný nebo experimentálně prokázaný účinek na protein nebo mRNA sestřih. Pokud poslední podmínka chybí, LR nejméně 350 vede nejvýše k PP1 Strong.

Produkční klasifikace je v `backend/classification_dag/`. Automatický výpočet
spouští `runtime.py`, jednotlivé rodiny pravidel jsou v podbalíku `nodes/` a společná
ENIGMA combinační politika je v `policy.py`. Soubor
`backend/modules/classifier.py` není produkční závislost. Dočasně slouží pouze
jako nezávislý oracle v regresních testech před jeho odstraněním. Společně s ním
se odstraní paritní testy a kompatibilní re-export
`backend/classification_dag/native.py`; produkční runtime tento soubor nepoužívá.

Transportní vrstva `backend/main.py` předává klasifikační požadavek službě
`backend/services/variant_classification_service.py`. Normalizaci, sestavení vstupu provider DAGu,
paralelní získání evidence, diagnostiku a fail-closed zpracování řídí
`backend/services/evidence_orchestration.py`. Veřejný Pydantic výsledek sestavuje
oddělený `backend/services/classification_presentation.py`. Prezentační služba
nesmí přidávat, potlačovat ani měnit sílu kritérií.

## 2. Zpracování vstupu

### 2.1 Normalizace vstupu a HGVS

Uživatel zadává gen a jednu variantu. Vstupní normalizační vrstva
`backend/modules/variant_input.py` přijímá:

- referenční transkriptovou notaci, například `c.303T>G`,
- notaci s accession prefixem, například `NM_007294.4:c.303T>G`,
- běžné kopírované varianty s oddělovací dvojtečkou nebo genovým prefixem,
  například `:c.303T>G`, `BRCA1:c.303T>G`, `BRCA1 c.303T>G` nebo
  `BRCA1 NM_007294.4:c.303T>G`,
- kombinovanou starší formu s `p.` následkem,
- genomickou variantu ve tvaru `chr17:43099813:C>T`, `17:43099813 C>T`
  nebo `17-43099813-C-T`.

Vstupní parser toleruje rozdíly, které nemění význam varianty. `NM`/`nm`,
`c.`/`C.`, `p.`/`P.`, nukleotidové symboly a operace `del`, `dup`, `ins`,
`delins` a `inv` mohou být zadány bez ohledu na velikost písmen. Mezery a
tabulátory jsou povoleny kolem dvojtečky a mezi `c.` a `p.` částí. Kombinovaná
notace může oddělit `c.` a `p.` část mezerou nebo lomítkem. Výstup se vždy vrací
v jednotném kanonickém formátu.

Pro genomickou variantu je povinná sestava `GRCh37` nebo `GRCh38`. Aplikace ji
neodhaduje. Genomický vstup se hledá v obousměrném indexu vytvořeném z
verzovaného coding SNV snapshotu a normalizovaného indel snapshotu. Jedna
jednoznačná shoda vrátí kanonickou `c.` notaci, třípísmennou `p.` notaci,
referenční transkript a zdroj normalizace. Nulová nebo víceznačná shoda ukončí
požadavek před klasifikací.

Po jednoznačném převodu genomické alely na `c.` notaci se `p.` následek vždy
znovu odvodí stejným sekvenčním HGVS enginem jako u přímého `c.` vstupu.
Protein uložený v souřadnicovém snapshotu slouží jen jako nezávislá kontrola.

#### Důležité omezení genomických vstupů

Transkriptové a genomické vstupy zatím nejsou zpracovány stejně obecně.
Přesný `c.` vstup pro podporovaný transkript může HGVS engine ověřit a převést
na `p.` přímo z lokální referenční sekvence, i když daná varianta není ve
snapshotu. Genomický vstup `g.` nebo `chr:pos:REF>ALT` se naproti tomu zatím
nejprve hledá v lokálním reverse indexu coding SNV a známých indelů. Varianta,
která v tomto indexu není, se odmítne, přestože by její následné `c. → p.`
zpracování engine zvládl.

Současný referenční balík proto deklaruje pouze schopnost
`reference_transcript_c_to_p`. cdot záznamy pro GRCh38 v balíku slouží HGVS
provideru, ale samy nezapínají obecný převod `g. → c.`. Balík neobsahuje lokální
genomovou FASTA ani úplnou produkční vrstvu pro obecné mapování GRCh37 a GRCh38.
Po úspěšném snapshotovém převodu na `c.` už obě vstupní cesty používají stejný
HGVS engine a stejný výpočet proteinového následku.

Budoucí sjednocení vyžaduje samostatnou vrstvu B s lokální genomovou sekvencí,
zarovnáními pro obě sestavy, validací gapped alignmentů, HGVS normalizací alely
a testy jednoznačnosti. Do té doby nelze rozsah podpory genomických vstupů
prezentovat jako ekvivalentní přímému `c.` vstupu.

Povolené transkripty jsou kontrolovány proti zvolenému genu. Aktuální registr
obsahuje `BRCA1: NM_007294.4` a `BRCA2: NM_000059.4`. Přidání dalšího genu
vyžaduje deklarovaný referenční transkript, protein a validovaný referenční balík,
nikoli ruční slovník jednotlivých variant.

Pokud uživatel uvede accession, verze je povinná a musí přesně odpovídat
schválenému transkriptu. `NM_007294` ani `NM_007294.3` se tiše nepřevádějí na
`NM_007294.4`. U holého `c.` vstupu je schválený transkript doplněn podle genu a
je vždy viditelně uveden ve výsledku. Tím zůstává implicitní transkript
auditovatelný a stará souřadnicová soustava nemůže projít jako jiná varianta bez
upozornění.

Samostatné API `POST /api/normalize` používá stejnou vrstvu jako klasifikační
endpoint. Vrací zadanou notaci, kanonickou `c.` a `p.` notaci, transkript,
sestavu a zdroj. Klasifikace proto normalizaci nemůže obejít.

#### Rozsah přijímaných zkrácených a alternativních zápisů

Formátová normalizace `c.` notace je obecná. Odstraňuje nadbytečné mezery,
sjednocuje `C.` na `c.`, převádí nukleotidy na velká písmena a operace `del`,
`dup`, `ins`, `delins` a `inv` na malá písmena. Přijme tedy například
`C. 5266 DUP c` jako `c.5266dupC`. Předpona transkriptu může být zapsána jako
`NM_...` i `nm_...`; mezery a tabulátory kolem dvojtečky jsou tolerovány.

U indelů engine parsuje zápis s uvedenou sekvencí i bez ní a ověří tvrzení proti
úplné sekvenci referenčního transkriptu. Například `c.5266dupC` a `c.5266dup`
kanonizuje na `c.5266dup`; `c.2102delA` a `c.2102del` kanonizuje na
`c.2102del`. Nesprávný suffix se odmítne jako neshoda reference. Toto chování
není omezeno na varianty uvedené v indelové mapě.

Podporované vstupy jsou tedy rozděleny takto:

| Vstupní odchylka | Zpracování |
| --- | --- |
| velikost písmen, mezery a tabulátory | obecná normalizace |
| úvodní `:` před `c.` | odstranění oddělovače |
| prefix genu, například `BRCA1:c.` nebo `BRCA1 c.` | prefix je autoritativní a formulář i backend přepnou na uvedený gen |
| podporovaný transkript `NM_007294.4` nebo `NM_000059.4` | určí BRCA1 nebo BRCA2 a musí mít správnou verzi |
| `c.` a `p.` oddělené mezerou nebo `/` | obecná normalizace |
| indel s uvedenou nebo vynechanou sekvencí | HGVS normalizace a ověření proti sekvenci |
| nesprávný suffix nebo neplatný indel | odmítnutí s důvodem před klasifikací |
| samotná `p.` notace bez `c.` nebo genomické souřadnice | nepodporováno, protože nemusí jednoznačně určit DNA variantu |

Parser nehledá první výskyt `c.` uvnitř libovolného textu. Pokud vstup výslovně
obsahuje `BRCA1`, `BRCA2` nebo podporovaný referenční transkript, tato informace
má přednost před aktuálním výběrem ve formuláři. Holá `c.` notace gen nemění.
Vzájemně rozporná kombinace, například `BRCA1 NM_000059.4:c.3703C>T`, neznámý
prefix nebo nesprávná verze transkriptu se odmítne s vysvětlením.

### 2.2 Kontrola referenční alely

Každá SNV prochází před klasifikací kontrolou referenční báze. Kontrola používá lokální verzovaná data, nikoli ruční seznam variant.

Zdroje kontroly:

1. coding SNV snapshot pro kódující pozice,
2. intronická souřadnicová mapa pro podporované intronické pozice.

Příklad:

| Vstup | Výsledek kontroly |
| --- | --- |
| `BRCA1 c.181T>G` | přijato, reference na `c.181` je `T` |
| `BRCA1 c.181A>C` | odmítnuto, reference na `c.181` není `A` |

Neshoda reference vrací HTTP 422 a klasifikace se nespustí. Pokud nainstalovaná data neumožňují referenci ověřit, aplikace postup ukončí místo tichého pokračování.

Implementace je v `backend/modules/reference_validation.py`.

### 2.3 Odvození a kontrola proteinového následku

Proteinový následek se za běhu počítá ze sekvence schváleného referenčního
transkriptu. Primárním enginem je `biocommons.hgvs==1.5.7`, providerem
`cdot==0.2.30` `JSONDataProvider` a zdrojem sekvencí je lokální
checksumovaný panelový balík. Kód je v `backend/modules/hgvs_engine.py`,
`backend/modules/hgvs_provider.py` a `backend/modules/panel_seqfetcher.py`.

Pro přesně popsané coding SNV, delece, duplikace, inserce a delins engine:

1. parsuje variantu na přesném accession s verzí,
2. ověří uvedenou referenční sekvenci,
3. použije HGVS pravidlo posunu k 3' konci, včetně variant přes hranice exonů,
4. mapuje změněnou transkriptovou sekvenci na schválený RefSeq protein,
5. vrátí kanonickou `c.` a třípísmennou `p.` notaci.

Tím se zpracují i varianty, které nejsou ve snapshotu. Například
`BRCA1 c.2102delA` i `c.2102del` vracejí `c.2102del
p.(Lys701SerfsTer2)`. Zadaný suffix se přijme pouze tehdy, pokud odpovídá
referenční sekvenci. Není k tomu ruční slovník variant.

Čistě intronická varianta, UTR varianta nebo exonová CNV s nejistými breakpointy
nemá z DNA zápisu jednoznačný proteinový produkt. V takovém případě je kanonický
výstup `p.?`. Veřejný výsledek vysvětluje, že následek nelze určit z DNA notace
samotné a že může být nutná RNA, transkriptová nebo breakpointová evidence.
`p.?` není zaměněno za chybu enginu ani za tvrzení, že protein zůstane beze změny.

Uživatel může `p.` notaci dodat, ale aplikace ji nepoužívá jako zdroj pravdy.
Porovná ji se sekvenčně odvozeným výsledkem. Rozpor vrátí HTTP 422 a klasifikace
se nespustí. Pokud je sekvenční následek `p.?`, konkrétní uživatelský proteinový
následek se bez validovaného RNA nebo breakpointového zdroje nepřijme.

Coding SNV a indel snapshoty jsou nezávislá regresní kontrola, zdroj souřadnic a
externích identifikátorů. Nejsou skrytým fallbackem pro výpočet proteinového
následku. Rozpor enginu s konkrétním známým následkem ve snapshotu ukončí
požadavek jako konflikt validovaných zdrojů. ENIGMA Table 9 zůstává zdrojem
funkční a RNA evidence, nikoli normalizačním slovníkem.

Uživatel může dodat také zkrácený legacy frameshift zápis, například
`p.(Cys1225fs)`. Zkrácený zápis je přijat pouze tehdy, když původní aminokyselina
a její pozice přesně souhlasí se sekvenčně odvozeným plným následkem, například
`p.(Cys1225SerfsTer10)`. ARIANE ve výsledku vždy vrátí plný kanonický zápis.
Zápis `p.(Arg1225fs)` ani `p.(Cys1226fs)` proto pro stejnou variantu neprojde.

ARIANE přijímá také běžný starší synonymní zápis, například `p.Val1653Val`
nebo `p.(Val1653Val)`, a převádí jej na současný kanonický zápis
`p.(Val1653=)`. Pravidlo je obecné pro všechny třípísmenné aminokyselinové
kódy a použije se pouze tehdy, když je aminokyselina před a za pozicí stejná.
Skutečná substituce, například `p.Val1653Ala`, se tímto pravidlem nezmění.

Příklad:

```text
BRCA1 c.303T>G
BRCA1 p.(Tyr101Ter)
variant_type = nonsense
```

Pokud uživatel zadá `p.` notaci, která není totožná ani povoleným zkráceným
frameshift ekvivalentem odvozeného následku pro referenční transkript, aplikace
vrátí HTTP 422. Rozpor se neřeší výběrem jedné z hodnot.

### 2.4 Odvození typu varianty

Typ varianty se odvozuje z kombinace `c.` a `p.` notace. Rozlišují se zejména:

- nonsense,
- frameshift,
- missense,
- synonymous,
- intronic,
- canonical splice site,
- in-frame deletion, insertion a delins,
- exon deletion a duplication,
- initiation codon,
- 5' UTR a 3' UTR.

Typ varianty určuje, které větve pravidel lze použít. Nonsense varianta bez proteinového následku se nesmí tiše považovat za missense.

## 3. Automaticky vyhodnocovaná kritéria

Síla kritéria se převádí na body:

| Směr | Síla | Body |
| --- | --- | ---: |
| patogenní | Very Strong | 8 |
| patogenní | Strong | 4 |
| patogenní | Moderate | 2 |
| patogenní | Supporting | 1 |
| benigní | Very Strong | -8 |
| benigní | Strong | -4 |
| benigní | Moderate | -2 |
| benigní | Supporting | -1 |

BA1 je samostatné benigní kritérium a vede přímo ke Class 1.

Ve veřejném výsledku jsou kritéria vždy řazena ve stejném pořadí jako v
ACMG/AMP a ENIGMA přehledu: `PVS1`, `PS`, `PM`, `PP`, `BA`, `BS`, `BP`.
Číslované kódy se uvnitř skupiny řadí vzestupně. Interní kvalifikátory, například
`PM5_PTC` nebo `BS1_Supporting`, zůstávají u svého základního kritéria. Stejné
řazení používá webová tabulka, batch výsledek, CSV export a výsledek manuální
revize.

### 3.1 BA1, BS1 a PM2 Supporting

Zdroj: lokální snapshoty gnomAD v2.1.1 exomes non-cancer a gnomAD v3.1.2 genomes non-cancer včetně pokrytí.

Pro BA1 a BS1 se používá výhradně nejvyšší non-cancer FAF95 v populacích AFR,
AMR, EAS, NFE a SAS. `popmax_AF` a běžné `AF` jsou pouze doplňující údaje.
Chybějící FAF95 se jimi nenahrazuje a uživatel dostane vysvětlení, proč BA1/BS1
nebylo možné vyhodnotit.

Těchto pět skupin odpovídá outbred non-founder populacím uvedeným v ENIGMA
v1.2 Appendix G: Non-Finnish European, African, Latino, East Asian a South
Asian. V2 skupiny ASJ a FIN a v3 skupiny AMI, ASJ a FIN jsou founder populační
kontext. MID a OTH nejsou součástí ENIGMA scoring setu, ale ARIANE je
neoznačuje automaticky za founder populace. AF, AC, AN a FAF95 všech těchto
vyloučených skupin se uchovávají a zobrazují v rozbalovací části `gnomAD
population frequency policy`. Nikdy nevstupují do maxima pro BA1/BS1.

Pro PM2 se přítomnost kontroluje pouze v pěti ENIGMA non-founder populacích.
Záznam přítomný pouze ve founder nebo jiné vyloučené populační skupině není
považován za přítomnost v outbred populaci. Jeho existence a hodnoty zůstávají
ve výsledném auditu. Pokud je záznam gnomAD filtrován, nelze z něj naopak
odvodit nepřítomnost a PM2 se nepoužije.

Před použitím BA1 nebo BS1 se kontroluje také filtr zdrojového záznamu a výjimka
pro známé patogenní founder varianty. Záznam, který neprošel filtrem gnomAD, se
pro BA1/BS1 nepoužije. Founder výjimka používá povinný verzovaný soubor
`backend/data/brca_pathogenic_founder_variants.json`. Obsahuje osm BRCA1/2
variant doložených v NCBI GeneReviews Table 6 nebo přímo záznamem ENIGMA
expert panelu pro BRCA1 c.181T>G. Každý záznam má zdroj, kontext, kanonickou
c. notaci a aliasy. Metadata obsahují checksum záznamů a checksumy zdrojových
stránek při sestavení. Poškozený nebo chybějící soubor nezpůsobí použití
BA1/BS1 bez kontroly. Kritérium se nepoužije a důvod se zobrazí uživateli.
Pokud varianta splní frekvenční práh a současně podléhá founder výjimce, ARIANE
vrátí splněný kód ve zvláštním poli `excluded_criteria`. Ve webovém výsledku je
odděleně uvedena dosažená síla, nula započtených bodů a konkrétní důvod
vyloučení. Taková položka se nikdy nepřenáší do `criteria`, součtu bodů ani
klasifikace. Například u `BRCA1 c.181T>G` se může zobrazit splněný práh
`BS1_Supporting`, ale BS1 se podle ENIGMA v1.2 nepoužije, protože jde o dobře
doloženou patogenní founder variantu.
Aktuální obsah lze znovu ověřit proti zdrojům příkazem
`python scripts/build_brca_founder_snapshot.py --check`. Builder vybírá sedm
řádků podle explicitního označení `Founder variant` v GeneReviews a samostatně
ověřuje BRCA1 c.181T>G proti záznamu ClinVar s klasifikací ENIGMA expert panelu.
ENIGMA v1.2 neposkytuje úplný strojově čitelný seznam všech founder variant na
světě. Snapshot je proto explicitní provozní seznam doložených výjimek, ne
odhad podle frekvence nebo původu populace. Novou výjimku lze přidat jen se
zdrojovým záznamem a novým checksumem.

| Kritérium | Podmínka |
| --- | --- |
| BA1 | frekvence nad 0,001 a průměrná hloubka alespoň 20 |
| BS1 Strong | frekvence nad 0,0001 a nejvýše 0,001, hloubka alespoň 20 |
| BS1 Supporting | frekvence nad 0,00002 a nejvýše 0,0001, hloubka alespoň 20 |
| PM2 Supporting | nepřítomnost v požadovaných non-cancer datasetech a průměrná hloubka alespoň 25 |

Jediné pozorování varianty v outbred populaci není podle ENIGMA informativní.
ARIANE v takovém případě nepoužije BA1, BS1 ani PM2. Pro BA1 nebo BS1 jsou
požadována alespoň dvě outbred pozorování v datasetu, který poskytl skórované
FAF95. Počet pozorování se bere z populačních AC, neodhaduje se z AF.

PM2 se nepoužívá pro malé indely do 50 bp. Vyloučení se kontroluje také přímo
podle operace v `c.` HGVS, nikoliv pouze podle odvozeného proteinového důsledku.
PTC vytvářející malý indel, například `BRCA1 c.5533_5534insG p.(Tyr1845Ter)`,
proto může vstoupit do PVS1/PM5 PTC větve, ale nesmí získat PM2.

Appendix G dovoluje PM2 Supporting pro větší inserce, delece a indely nad 50 bp,
pokud jsou nepřítomné ve vhodném datasetu, prošly kontrolou kvality a byly
porovnány zahrnuté exony. Exonové CNV s neurčitými breakpointy se proto
nehledají jako přesné VCF ID. Obecná větev nejprve určí exon z Table 4, načte
jeho reprodukovatelně odvozený GRCh37 interval a v úplném verzovaném gnomAD-SV
datasetu hledá deleci zahrnující celý kódující interval exonu. Manifest
neobsahuje žádné konkrétní varianty ani předem přidělená kritéria.
Chybějící záznam bez odpovídajícího detekčního rozsahu nebo QC není důkazem
nepřítomnosti. Fixture ani neúplná cache nemůže vytvořit frekvenční kritérium.

ENIGMA v1.2 požaduje průměrnou hloubku v oblasti varianty, ale neurčuje počet
okolních bází. ARIANE proto nepoužívá vlastní pevné okno ±N bp. Oblast je
definována jako úplný genomový rozsah referenční alely varianty. U SNV jde o
jednu bázi, u vícenukleotidové substituce o všechny báze `REF`. Hloubka je
aritmetický průměr per-position hodnot. Chybí-li jediná pozice rozsahu, pokrytí
není prokázáno a PM2 se nepoužije. Ve výsledných datech jsou uloženy hranice
rozsahu, očekávaný a dostupný počet pozic a použité klíče coverage snapshotu.

### 3.2 PVS1 a PM5 PTC

Zdroj: ENIGMA Specifications Table 4 v1.2.

PVS1 se vyhodnocuje pro:

- nonsense varianty,
- frameshift varianty,
- vybrané splice-site varianty,
- exonové delece,
- exonové duplikace.

Table 4 obsahuje pravidla pro jednotlivé exony, kritické C-terminální hranice, splice varianty a exonové přestavby. Výsledná síla může být Very Strong, Strong, Moderate, Supporting nebo N/A.

U splice variant se PVS1 nepřidává pouze podle vzdálenosti od exonu. Varianta musí mít odpovídající pravidlo v Table 4. Větve závislé na RNA se automaticky použijí jen při přesné shodě s úplnou ENIGMA Supplementary Table 2 a při jednoznačném průchodu pravidly Appendix E. Ostatní RNA větve zůstávají bez bodů a jsou označeny k odborné revizi.

U iniciačního kodonu se automatické PVS1 nepoužívá. Aplikace vytvoří doporučení pro strukturovanou revizi podle iniciačního flowchartu.

PVS1 a PM5 PTC se vždy přebírají ze stejného vybraného řádku Table 4. U
frameshiftu určuje Appendix D sílu PM5 podle exonu, ve kterém leží nukleotidová
změna, i když se nový terminační kodon nachází v pozdějším exonu. Proteinový
následek a pozice nového stop kodonu se nadále přesně počítají a zobrazují, ale
nepoužívají se k nezávislému přepnutí PM5 do jiného řádku Table 4.

Poznámka ke zdroji: stručný Readme v souboru Table 4 popisuje exon terminačního
kodonu, zatímco podrobný Appendix D na straně 27 výslovně určuje exon
nukleotidové změny a uvádí frameshift se stop kodonem v následujícím exonu.
Runtime se řídí podrobným postupem a příkladem z Appendix D. Původní Readme je
ve zdrojovém snapshotu zachován beze změny pro audit.

V koncové hraniční větvi se řádek vybírá podle začátku změněné nebo zkrácené
kritické proteinové sekvence v souladu s PVS1 flowchartem. Například `BRCA1
c.5556_5560del p.(Tyr1853AspfsTer25)` používá řádek PVS1 Very Strong a PM5 PTC
Strong. Výsledkem je 12 bodů a třída 5.

#### 3.2.1 PVS1 z RNA evidence

Automatické `PVS1_RNA` nevychází z ručního seznamu variant. Klasifikátor používá
úplnou oficiální ENIGMA Supplementary Table 2 a obecnou rozhodovací větev z
Appendix E Table 9. Kritérium se přidělí pouze tehdy, když jsou současně splněny
všechny následující podmínky:

- přesná normalizovaná varianta je obsažena v Supplementary Table 2;
- ENIGMA ji zařadila do kategorie pacientské mRNA bez alelově specifické
  kvantifikace s aberantními transkripty odpovídajícími ztrátě funkce;
- výsledek jednoznačně popisuje deleci jednoho celého exonu;
- uvedený exon se jednoznačně mapuje na deleční řádek Table 4;
- Table 4 pro tento transkriptový důsledek uvádí aplikovatelnou základní PVS1
  sílu.

Kvalitativní větev pacientské mRNA bez alelově specifické kvantifikace snižuje
základní sílu podle Appendix E: Very Strong na Strong, Strong na Moderate a
Moderate na Supporting. Komplexní nebo částečné transkriptové následky se
neodhadují. Stejně tak se bez uloženého procenta funkčního transkriptu
nepoužijí kvantitativní větve Appendixu.

Při přijetí `PVS1_RNA` se podle Figure 1B odstraní slabší predikční evidence pro
stejný splice mechanismus, například PP3, BP4, BP7, BP1 nebo predikční PS1.
Proteinová funkční evidence PS3/BS3 se automaticky nemaže, protože může
popisovat jiný mechanismus; případná kombinace se označí k odborné kontrole.

Například `BRCA1 c.4185G>A` má v oficiální Supplementary Table 2 pacientskou
mRNA s delecí exonu 12. Table 4 uvádí pro deleci BRCA1 E11(12) základní PVS1.
Kvalitativní RNA větev proto vrátí `PVS1 Strong (RNA)`.

### 3.3 PS3 a BS3

Zdroj: ENIGMA Specifications Table 9 v1.2.

Vyhledávání používá přesný klíč `gene:c_notation`. Automaticky se použijí pouze řádky s přiřazeným PS3 nebo BS3 a podporovanou silou.

Table 9 obsahuje také řádky, ve kterých PS3 ani BS3 nebylo splněno. Tyto řádky zůstávají součástí lossless snapshotu, ale nevytvářejí kritérium.

Sloupec `Splicing Prediction` v Table 9 zaznamenává nejvyšší ze čtyř delta
skóre SpliceAI s oknem 10 kb, které ENIGMA použila jako kontext při posouzení
funkčních testů pro PS3 a BS3. Tato hodnota nepřepisuje SpliceAI výsledek
konfigurovaného zdroje pro Figure 1A, tedy pro PP3, BP4, BP7 a BP1, ani pro
splice podmínku proteinového PS1 u hodnocené varianty.

Pokud se obě hodnoty liší, ARIANE zobrazí obě hodnoty a jejich původ. Jestliže
leží v různých ENIGMA pásmech `<= 0,1`, `> 0,1 a < 0,2` nebo `>= 0,2`, výstup
vyžaduje odbornou kontrolu provenance SpliceAI. Table 9 se přesto nadále použije
pro své explicitní doporučení PS3 nebo BS3 a pro publikovanou RNA informaci.
Pokud konfigurovaný SpliceAI výsledek chybí, hodnota z Table 9 jej nenahrazuje a
závislá automatická kritéria zůstávají nedostupná.

Funkční evidence se do automatické klasifikace nepřidává z tutorialu ani z
variantově specifického lokálního záznamu. Varianta získá PS3 nebo BS3 pouze
tehdy, když přesný řádek úplné ENIGMA Table 9 obsahuje aplikovatelné doporučení.
Jinak se kritérium nepřidělí a případná další publikovaná evidence patří do
odborné revize podle Figure 1C a Appendix E.

PS3 a BS3 nemají dva oddělené lookup moduly. Obě kritéria jsou opačné výsledky
stejného kalibrovaného funkčního hodnocení a proto je společně vrací
`backend/modules/table9.py`; v klasifikačním DAG je zpracovává uzel
`rule.functional.table9`. Samostatná BP7 RNA větev tento výsledek pouze čte,
nevytváří náhradní BS3.

### 3.4 PP4 a BP5

Zdroj automatických kritérií: verzovaný snapshot variantově specifických combined clinical LR odvozený ze dvou veřejných ENIGMA zdrojů. Multifaktoriální klinické LR pocházejí z UCSC ENIGMA `BRCAmfa` tracku. Case-control LR pocházejí ze Supplementary Data 5 studie Zanti et al. 2025, která vznikla v ENIGMA Analytical Working Group. Supplementary Table 7 ani její posterior probability se pro PP4/BP5 nepoužívají.

| Kritérium | Supporting | Moderate | Strong | Very Strong |
| --- | ---: | ---: | ---: | ---: |
| PP4 | LR >= 2,08 | LR >= 4,3 | LR >= 18,7 | LR >= 350 |
| BP5 | LR <= 0,48 | LR <= 0,23 | LR <= 0,05 | LR <= 0,00285 |

PP4 a BP5 se automaticky vyhodnocují z lokálního verzovaného snapshotu variantově specifických klinických LR. Manuální revize zůstává dostupná pro varianty nebo zdroje, které snapshot neobsahuje. Reviewer zadá variantově specifickou klinickou hodnotu, její škálu, citaci zdroje a souhrn zahrnutých klinických dat včetně kontroly jejich nezávislosti. Podporované škály jsou běžný LR, `log10(LR)` a ACMG evidence points. ARIANE určí sílu výhradně podle ekvivalentních prahů. Jedna publikace stačí, pokud poskytuje metodicky přijatelný variantově specifický klinický LR. Není nutné kombinovat více publikací. Sílu PP4 nelze ručně přepsat a neúplný záznam nelze aplikovat.

ENIGMA v1.2 požaduje pro PP4/BP5 combined LR klinických dat. Historická
multifaktoriální posterior probability ani IARC třída se nepoužívají přímo jako
PP4/BP5 a nepřevádějí se pomocí pevného obecného prioru. Tutorialový nebo starší
výsledek založený na posterior probability se proto může lišit od aplikace
současných v1.2 LR prahů.

Automatický snapshot je uložen v souborech:

- `data/precomputed/brca_pp4_clinical_lr_snapshot.index.json`,
- `data/precomputed/brca_pp4_clinical_lr_snapshot.metadata.json`.

Builder `scripts/build_pp4_clinical_lr_snapshot.py` používá všechny čtyři variantově specifické komponenty, ze kterých oficiální UCSC ENIGMA `BRCAmfa` track verze 1.1.0 přepočítává combined LR: Easton et al. 2007, PMID 17924331; Parsons et al. 2019, PMID 31131967; Li et al. 2020, PMID 31853058; a Caputo et al. 2021, PMID 34597585. K nim přidává case-control LR ze studie Zanti et al. 2025, PMID 40413188, DOI `10.1038/s41467-025-59979-6`.

Ze Zanti Supplementary Data 5 se přijímá publikovaný finální soubor 1 710 variant: 681 BRCA1 a 1 029 BRCA2. Builder vyžaduje oblast CDS ±5 bp, non-founder FAF nejvýše 0,001, alespoň tři nositele v kombinovaných datech a u BRCA2 evidence alespoň ze dvou datasetů. Dvě další BRCA2 varianty mají v publikaci LR `N/A`, protože výpočet nekonvergoval. Nejsou bodovány a zůstávají výslovně uvedené v metadatech. Snapshot používá přímo klinické LR. Posteriorní pravděpodobnost se nepřevádí pomocí obecného prioru.

Před uložením snapshotu projde každá zdrojová `c.` notace stejným lokálním
`biocommons.hgvs` enginem jako uživatelský vstup. Uvedená sekvence delece nebo
duplikace se ověří proti checksumovanému referenčnímu transkriptu a uloží se
kanonická notace i ověřené zdrojové aliasy. Známé indely se navíc křížově
kontrolují proti normalizovanému indelovému snapshotu. Nejde o runtime fallback
ani o ruční slovník. Nevalidní notace a konflikt normalizovaných zdrojů se
nezařadí a důvod zůstane v metadatech.

Pokud více řádků po normalizaci popisuje stejnou alelu, builder spojí její přijaté komponenty a přepočítá jediný combined LR. Kritérium jednotlivé komponenty se samostatně neboduje. Z výsledného LR vznikne právě jedno PP4 nebo BP5. Stejný `source_id` ani stejná skupina nezávislosti se nesmí započítat dvakrát. Duplicita je fatální chyba buildu. Rozdílné hodnoty pod stejným zdrojem jsou konflikt a daná alela se automaticky nepoužije.

Metadata obsahují zdrojový manifest, URL a checksum obou zdrojových datasetů, checksum indexu, verzi pravidel,
provenance HGVS enginu a referenčního balíku, checksum závislého indelového
snapshotu, počty záznamů a seznam konfliktů. Chybějící metadata, nesprávný
checksum indexu nebo zdrojového manifestu, nesprávný počet záznamů nebo nejednoznačný alias zastaví spuštění aplikace. Pro `BRCA1 c.5266dup` a alias `c.5266dupC` je combined LR z Li et al. 2020 a Zanti et al. 2025 `1,36181 × 10^90`, což odpovídá PP4 Very Strong. Pro
`BRCA2 c.9891_9894dup` a zdrojový zápis `c.9891_9894dupATTT` obsahuje LR
`0,41018` ze studie Li et al. 2020, což odpovídá BP5 Supporting.

Pro `BRCA1 c.509G>A` se multifaktoriální LR `6,1764` násobí case-control LR `0,00639025`. Výsledný combined LR je `0,0394687`, proto se aplikuje jediné kritérium BP5 Strong. Dílčí PP4 Moderate se samostatně neaplikuje.

Síla BP5 a použití této síly v klasifikační kombinaci jsou dvě oddělená
rozhodnutí. ENIGMA Table 3 dovoluje přiřadit Likely Benign z jediného Strong
benigního kódu pouze tehdy, když do něj přispělo více typů evidence. ARIANE
proto u každého BP5 Strong zachovává počet LR příspěvků a seznam klinických
typů evidence. Jediné BP5 Strong může samo vést ke Class 2 pouze při nejméně
dvou zaznamenaných LR příspěvcích ze dvou klinických typů evidence. Jeden LR
zůstává platným BP5 Strong, ale bez další benigní evidence sám nestačí na
Class 2. Aktuální snapshot obsahuje 378 BP5 Strong záznamů. Z nich 218 má
nejméně dva zaznamenané LR příspěvky a 160 má pouze jeden LR příspěvek.

Aktuální snapshot obsahuje 5 147 jednoznačných variantových záznamů. Zanti case-control komponentu obsahuje 1 710 záznamů.

Zdrojový manifest eviduje všechny studie vyjmenované v ENIGMA Appendix B a navíc Caputo 2021 a Zanti 2025, jejichž variantově specifická data používá automatický výpočet:

| Zdroj | PMID |
| --- | --- |
| Goldgar et al. | 15290653 |
| Thompson et al. | 12900794 |
| Easton et al. | 17924331 |
| Spurdle et al. | 25857409 |
| de la Hoya et al. | 27008870 |
| Parsons et al. | 31131967 |
| Li et al. | 31853058 |
| Caputo et al. | 34597585 |
| Zanti et al. | 40413188 |

Status zdroje se zpracovává fail-closed:

- `ENIGMA recognised source`: musí být vybrán jeden z uvedených PMID. Po splnění ostatních požadavků může PP4 vstoupit do amended klasifikace.
- `Other reviewed source`: vyžaduje citaci, jméno reviewera a metodické zdůvodnění kompatibility s ENIGMA PP4. Po splnění ostatních požadavků může PP4 vstoupit do amended klasifikace.
- `Unreviewed source`: hodnota a zdroj se zachovají v auditním záznamu, ale PP4 se neaplikuje a nepřidají se body.

Za primární zdroj evidence se považuje publikace nebo verzovaný dataset. Sekundární databáze, například CANVarUK, může sloužit k nalezení a zobrazení hodnoty, ale nenahrazuje citaci primárního zdroje. Pokud amended výsledek kombinuje PP1 nebo PS4 s automatickým či manuálním PP4/BP5, backend vyžaduje explicitní potvrzení nezávislosti pozorování a textové zdůvodnění. Bez nich požadavek odmítne. Samotné zaškrtnutí nenahrazuje kontrolu zdrojových kohort a klinických LR komponent.

### 3.5 PS1 na proteinové úrovni

Proteinové PS1 používá oficiální ENIGMA Supplementary Table 7 jako referenční
P/LP klasifikační dataset. Z jeho 146 P/LP variant je 60 normalizovaných
missense referencí zařazeno do proteinového PS1 registru. Každý záznam má
samostatný stav `eligible`, `excluded` nebo `review_required`, protože samotná
P/LP klasifikace ještě nenahrazuje PS1-specifickou splice kontrolu.

PS1 vyžaduje:

- missense variantu,
- stejnou normalizovanou missense substituci jako známá P/LP varianta,
- jinou nukleotidovou změnu,
- P/LP klasifikaci reference v ENIGMA ST7 v1.2, oficiální ENIGMA/ClinGen VCEP
  assertion nebo úplnou lokální reklasifikaci podle uvedené verze ENIGMA VCEP,
- SpliceAI nejvýše 0,1 u reference i hodnocené varianty,
- žádný potvrzený škodlivý splice efekt u obou variant po kontrole definovaných
  verzovaných ENIGMA zdrojů.

Patogenní referenční varianta dává PS1 Strong. Likely Pathogenic referenční varianta dává PS1 Moderate.

Proteinové PS1 používá pro hodnocenou variantu SpliceAI výsledek z
konfigurovaného ENIGMA-kompatibilního zdroje. Hodnota `spliceai_prediction` z
Table 9 je auditní kontext funkčního přezkumu PS3/BS3 a tento výsledek
nenahrazuje. Pokud je konfigurovaný výsledek nedostupný, proteinové PS1 se
automaticky nepřidělí. Rozdíl proti Table 9 se zobrazí a rozdíl mezi ENIGMA
predikčními pásmy vyžaduje odbornou kontrolu provenance.

Automatické body lze přidat pouze ze záznamu se stavem `eligible` v
`backend/data/ps1_protein_reference_registry.json`. Registr nyní obsahuje 60
ST7 missense referencí: 37 `eligible` a 23 `excluded`. Každý záznam obsahuje
původ klasifikace, splice stav reference, SpliceAI provenance, kontrolované
zdroje, datum, checksum podkladu a známé PS1 závislosti. Validátor odmítá přímou
i delší známou kruhovou závislost.

Záznam `review_required` se zobrazí jako předvyplněný kandidát pro
strukturovanou manuální revizi `PS1_PROTEIN`, bez bodů. Záznam `excluded` se
zobrazí s důvodem vyloučení a nelze jej manuálně potvrdit jako proteinové PS1.
Potvrzený nebo predikovaný splice efekt proteinové PS1 vylučuje. Konfliktní
nebo neúplná evidence vede k revizi.

Stav `none_identified` neznamená, že splice efekt neexistuje. Znamená pouze, že
nebyl nalezen při kontrole přesně uvedených verzí ENIGMA Table 9 a
Supplementary Table 2.

PS1 pro stejný splice efekt se automaticky neboduje. Aplikace pouze označí kandidáta pro strukturovanou manuální revizi.

Manuální formulář pro proteinové PS1 přijímá referenční c. HGVS. Backend
referenci normalizuje proti nakonfigurovanému transkriptu, odvodí a ověří její
p. následek, porovná jej s hodnocenou variantou a získá SpliceAI pro obě
varianty ze stejného profilu. Současně ověří přesnou referenci v ClinVar a
ClinGen ERepo. Frontend pouze zobrazí vrácená fakta.

Hvězdičky ClinVar popisují review status a žádný jejich počet nekvalifikuje
referenci pro PS1. Samostatně identifikovaná assertion příslušného ENIGMA/ClinGen
VCEP může být způsobilá nezávisle na aggregate star rating. Běžný ClinVar
aggregate závěr nepředvyplňuje klasifikaci reference, její ověření, klasifikační
zdroj ani evidenční reference. Jinak by se závěr klasifikace znovu použil jako
její vlastní podklad. Ani u skutečné VCEP assertion se nepřeskakuje kontrola
definovaných RNA/splice zdrojů a přímé reciproční PS1 závislosti. Nedostupné p.
odvození nebo SpliceAI zůstane nedostupné a PS1 se nepřidá.

### 3.6 PP3

PP3 Supporting má dvě automatické větve.

SpliceAI větev:

- SpliceAI alespoň 0,2,
- pouze pro povolené typy, například synonymous, missense, in-frame a intronic,
- nepoužívá se jako obecné PVS1 pro nonsense, frameshift, exonové CNV nebo canonical splice-site varianty.

BayesDel_noAF větev:

- missense nebo in-frame varianta v klinicky významné funkční doméně,
- BRCA1: BayesDel_noAF alespoň 0,28,
- BRCA2: BayesDel_noAF alespoň 0,30.

BayesDel se vyhodnocuje pouze po skutečném výsledku SpliceAI a průchodu
příslušnou větví Figure 1A. Chybějící SpliceAI není pásmo `no impact` ani
`not informative`. Při nedostupném SpliceAI se výsledek Figure 1A označí jako
nedostupný a nepoužije se PP3, BP4, BP1 ani BP7. Skutečně naměřené SpliceAI
větší než 0,1 a menší než 0,2 zůstává oficiálním pásmem `not informative`.
U missense nebo in-frame varianty může tato větev podle Figure 1A pokračovat
přes funkční doménu k BayesDel.

PP3 se nepřičítá současně s aplikovaným PVS1.

### 3.7 BP4

BP4 Supporting vyžaduje potvrzený nízký SpliceAI. Chybějící skóre se nepovažuje za nulové.

Pro missense a in-frame variantu uvnitř funkční domény platí:

- SpliceAI nejvýše 0,1,
- BRCA1 BayesDel_noAF nejvýše 0,15,
- BRCA2 BayesDel_noAF nejvýše 0,18.

Pro synonymous variantu uvnitř domény a pro podporovanou intronickou variantu se používá potvrzené SpliceAI nejvýše 0,1 podle příslušné větve pravidla.

### 3.8 BP7

BP7 Supporting se používá společně s BP4.

U synonymous variant uvnitř funkční domény vyžaduje aplikované BP4 a SpliceAI nejvýše 0,1. U intronických variant se navíc kontroluje, že pozice neleží v konzervovaném donorovém nebo akceptorovém motivu.

Synonymous varianta mimo funkční doménu je řešena přes BP1, nikoli přidáním BP7.

RNA větev BP7 Strong vyžaduje manuální strukturovanou revizi.

### 3.9 BP1

BP1 Strong se používá pro missense, synonymous a in-frame varianty mimo klinicky významnou funkční doménu, pokud je SpliceAI potvrzeně nejvýše 0,1.

Použité domény:

| Gen | Doména | Aminokyselinový rozsah |
| --- | --- | --- |
| BRCA1 | RING | 2 až 101 |
| BRCA1 | coiled-coil | 1391 až 1424 |
| BRCA1 | BRCT | 1650 až 1857 |
| BRCA2 | PALB2 binding | 10 až 40 |
| BRCA2 | DBD | 2481 až 3186 |

Chybějící SpliceAI znamená, že BP1 nelze použít.

Příklad: `BRCA1 c.3247A>G p.(Met1083Val)` leží mimo uvedené domény. Broad
SpliceAI s profilem Appendix J pro referenční transkript `ENST00000357654.9`,
který odpovídá `NM_007294.4`, vrací maximum raw delta `0,01`. Varianta proto
splňuje BP1 Strong.

## 4. Kritéria vyžadující manuální revizi

Automatická Module 1 klasifikace nepřidává PS4, PM3, PP1, BS2 a BS4. PP4 a BP5 přidává pouze při přesné shodě s validovaným lokálním snapshotem klinických LR. PVS1 RNA může přidat z přesného oficiálního ST2 záznamu a obecného rozhodovacího pravidla popsaného výše. Ostatní uvedené kódy závisejí na datech, která nelze bezpečně odvodit pouze z HGVS varianty.

Strukturovaná manuální část dále podporuje:

- PP4 z variantově specifického combined clinical LR,
- PVS1 RNA pro další publikovanou nebo komplexní RNA evidenci, kterou nelze
  jednoznačně vyhodnotit z úplné ST2 a Table 4,
- BP7 RNA,
- PVS1 pro iniciační kodon,
- PS1 pro stejný splice efekt.

Manuálně doplněná kritéria vytvářejí oddělený amended working result. Původní automatická Module 1 klasifikace zůstává zachována.

## 5. Postup klasifikace

Po normalizaci vstupu vytvoří backend `ClassificationRequest` a předá jej
produkčnímu grafu `4.0.0-gene-policy-provider-dag`. `main.py` klasifikační zdroje přímo
nevolá. Provider uzly grafu získají:

- souřadnice GRCh37 a GRCh38,
- SpliceAI pro hodnocenou variantu a potřebné kandidátní PS1 reference,
- BayesDel_noAF a informační AlphaMissense,
- gnomAD,
- ENIGMA Table 9,
- combined clinical LR,
- exonovou CNV evidenci,
- kurátorovanou informaci o důležitém reziduu,
- podklady proteinového PS1.

Nezávislé providery běží paralelně. Každý vrací typovaný `EvidenceItem` se stavem,
důvodem a provenance. `EvidenceBundle` vznikne až po dokončení provider vrstvy.
Nedostupná hodnota zůstává `UNAVAILABLE` a nepřevádí se na nulu ani na negativní
výsledek pravidla. ClinVar a ClinGen ERepo jsou pouze externí porovnání a nejsou
vstupem klasifikace.

Kritéria se po sestavení evidence vyhodnocují v tomto logickém pořadí:

1. kontrola vstupu, reference a proteinového následku,
2. souřadnice GRCh37 a GRCh38, pokud má varianta jednoznačně definovatelnou
   genomovou alelu,
3. gnomAD a BA1,
4. Table 9 pro PS3 a BS3,
5. Table 4 pro PVS1 a PM5 PTC,
6. lokální snapshot klinických LR pro PP4 nebo BP5,
7. proteinové PS1,
8. PP3 a BP4,
9. BP7,
10. BP1,
11. klasifikační kombinace.

Exonové delece a duplikace s neurčitými breakpointy nepokračují do
souřadnicového resolveru, SpliceAI ani BayesDel. Tyto služby vyžadují konkrétní
genomovou alelu a jejich použití na intervalový exonový zápis by nebylo
interpretovatelné. PVS1 se pro takovou variantu vyhodnotí přímo z ENIGMA Table 4.

BA1 ukončí klasifikaci jako Class 1.

Pokud jsou všechna aplikovaná kritéria pouze v jednom směru, používají se kombinace ENIGMA VCEP v1.2 Table 3. Samotný součet bodů v tomto případě neurčuje třídu. Tavtigian 2020 je bodový systém a pro evidenci pouze jedním směrem se v ARIANE nepoužívá.

Jedno Very Strong benigní kritérium postačuje podle klasifikačního textu ENIGMA
pro Likely Benign. Jedno Strong benigní kritérium postačuje pouze tehdy, pokud
konkrétní záznam kritéria auditovatelně dokládá více příspěvků evidence. Tato
vlastnost není odvozena z názvu kódu:

- BP1 Strong je způsobilé po úplném průchodu Figure 1A, protože kombinuje typ
  varianty, polohu vůči funkční doméně a negativní splice predikci.
- BS4 Strong je způsobilé pouze při nejméně dvou doložených a nezávisle
  identifikovaných segregačních LR příspěvcích.
- BP5 Strong je způsobilé pouze při nejméně dvou zaznamenaných LR příspěvcích
  ze dvou klinických typů evidence.
- Chybějící nebo neúplná provenance znamená, že samotné Strong kritérium
  nestačí na Likely Benign. Síla kritéria ani jeho body se tím nemění.

Pokud se Strong kritérium kombinuje s další benigní evidencí podle Table 3,
platí běžné kombinace Table 3 bez ohledu na tuto výjimku pro jediné Strong.
Při mixed evidence se nadále použije druhý klasifikační postup a bodový systém.

Příklad: PVS1 Very Strong bez dalšího kritéria zůstává Class 3, protože nesplňuje kombinaci pro Likely Pathogenic.

Síly patogenní evidence se při použití Table 3 vyhodnocují monotónně. Kritérium
Very Strong může splnit požadavek na další Strong kritérium. Dvě nezávislá Very
Strong kritéria proto vedou k Class 5. Silnější evidence nesmí snížit třídu, které
by bylo dosaženo se stejným kritériem na úrovni Strong. Regresním příkladem je
BRCA1 c.4676-1G>A s PVS1 Very Strong a PP4 Very Strong.

Stejná monotónní zásada platí pro benigní směr. Zesílení nebo přidání benigní
evidence nesmí posunout výsledek k vyšší třídě. Testovací matice vyčerpávajícím
způsobem kontroluje přidávání a postupné zesilování kritérií v obou směrech.

Pokud jsou současně přítomna patogenní i benigní kritéria, nastává druhý klasifikační postup ENIGMA. V tomto případě se používá bodový systém Tavtigian 2020:

| Součet | Třída |
| ---: | --- |
| 10 a více | Class 5, Pathogenic |
| 6 až 9 | Class 4, Likely Pathogenic |
| -1 až 5 | Class 3, VUS |
| -6 až -2 | Class 2, Likely Benign |
| méně než -6 | Class 1, Benign |

Výsledek s protichůdnými směry zachovává vypočtenou ENIGMA třídu a obsahuje barevný pruh `Mixed evidence`. Pruh uvádí, že byla použita ENIGMA bodová kombinace a že je nutná expertní revize. Odkazuje přímo na verzovaný dokument [ENIGMA Specifications v1.2](https://cspec.genome.network/cspec/File/id/11e62fec-23b0-4a3e-b2df-751855301746/data), část `Classification Methods`, druhý postup. Rozbalovací technický detail zvlášť ukazuje součet patogenních bodů, benigních bodů a celkový výsledek.

ENIGMA v úvodní poznámce Specifications v1.2 upozorňuje, že pravidla nemusí
spolehlivě odlišit varianty se středním rizikem nebo sníženou penetrancí. Rozpor
mezi více typy evidence má vést k dalšímu zkoumání snížené penetrance nebo
částečného účinku na funkci či splicing. ARIANE proto při souběhu funkčního PS3
a benigního klinického LR BP5 zobrazí upozornění k expertní revizi. Samotný
konflikt není důkazem snížené penetrance a nemění kritéria ani výpočet třídy.

## 6. Oficiální ENIGMA tabulky

Oficiální zdroje jsou genově oddělené záznamy ClinGen CSpec pro ENIGMA BRCA1/2 VCEP v1.2.0, vydání 2025-01-09:

- BRCA1: [ClinGen CSpec GN092](https://cspec.genome.network/cspec/ui/svi/doc/GN092)
- BRCA2: [ClinGen CSpec GN097](https://cspec.genome.network/cspec/ui/svi/doc/GN097)

### 6.0 Kompletní veřejný prohlížeč tabulek

Oficiální balík obsahuje tři samostatné číselné řady, které se nesmějí zaměňovat:

- 9 Specification Tables;
- 17 Appendix Tables;
- 16 Supplementary Tables.

Všech 42 tabulek je reprodukovatelně převedeno do prezentačního souboru
`backend/data/enigma_reference_tables.json`. Generátor
`scripts/build_enigma_reference_tables.py` čte přímo checksumem ověřené DOCX a
XLSX soubory v `docs/enigma/v1.2/source`. U každé tabulky ukládá řadu, číslo,
název, části nebo listy, zdrojové řádky, vzorce buněk, identifikátor zdroje a
checksum. Obsah variant ani rozhodovací hodnoty se nepřepisují ručně.

Tento soubor slouží veřejnému prohlížeči a přesně vymezenému lookupu klinických
anotací z Appendix Table 11. Klasifikační kritéria dál používají menší účelové
runtime datasety s vlastními validačními kontrolami, například Table 4, Table 9,
ST2 a ST7. Rozhraní u každé tabulky rozlišuje přímý runtime lookup,
implementovanou definici pravidla, expertní revizi, kandidátní registr a
kalibrační nebo referenční podklad.

Prohlížeč tabulek řadí obsah podle role, nikoliv pouze podle dokumentu:

- `Used by ARIANE`, 15 tabulek: Specification Tables 1, 2, 3, 4, 7 a 9;
  Appendix Tables 3, 4, 9, 11, 14, 15 a 16; Supplementary Tables 2 a 7;
- `Expert review`, 4 tabulky;
- `Supporting and calibration`, 24 tabulek.

Označení `Used by ARIANE` zahrnuje přímé runtime lookupy, implementované definice
pravidel a schválené kandidátní registry. Neznamená, že se všech 15 tabulek při
každé klasifikaci načítá jako jeden obecný sešit. Každé pravidlo používá účelový
validovaný dataset nebo explicitní implementaci odpovídající příslušné tabulce.

Strukturovaný pohled zachovává hodnoty a vzorce, ale může zjednodušit sloučené
buňky a původní vzhled. Autoritativní vizuální podoba zůstává v odkazovaném
oficiálním dokumentu.

### 6.1 Table 4

Runtime soubor: `backend/data/enigma_table4.json`

Generátor: `scripts/build_enigma_table4_snapshot.py`

Obsah:

- 493 zdrojových řádků,
- 20 zdrojových sloupců,
- indexy exonových rozsahů,
- PTC pravidla,
- splice pravidla,
- pravidla exonových delecí,
- pravidla exonových duplikací.

### 6.2 Table 9

Runtime soubor: `backend/data/enigma_table9.json`

Generátor: `scripts/build_enigma_table9_snapshot.py`

Obsah:

- 4 731 řádků,
- 14 zdrojových sloupců,
- PS3 a BS3 přiřazení,
- publikované splice výsledky,
- revidované SpliceAI hodnoty,
- 437 revidovaných řádků bez aplikovaného PS3 nebo BS3.

### 6.3 Supplementary Table 7

Runtime soubor: `backend/data/st7_reference_set.json`

Generátor: `scripts/build_enigma_st7_snapshot.py`

Obsah:

- 773 variant,
- 28 zdrojových sloupců,
- prior a posterior probability,
- IARC třída,
- populační a referenční údaje.

P/LP missense záznamy ST7 se používají jako referenční klasifikační základ
proteinového PS1. O automatickém použití rozhoduje až explicitní stav v
proteinovém PS1 registru po kontrole SpliceAI, Table 9 a úplné ST2.

### 6.3.1 Supplementary Table 2, PVS1 RNA a proteinové PS1

Úplný runtime soubor: `backend/data/enigma_st2_splice_evidence.json`

Generátor: `scripts/build_enigma_st2_splice_evidence_snapshot.py`

Obsahuje všech 220 variant a všech 11 zdrojových sloupců listu
`ST2 splicing dataset codes`, číslo zdrojového řádku a checksum oficiálního
Excelu. Používá se ke kontrole známé RNA/splice evidence pro proteinové PS1 a
jako oficiální vstup obecné automatické větve PVS1 RNA popsané v části 3.2.1.
ARIANE nemá aktivní registr referencí pro splice PS1. Úplná ST2 se používá
jen jako definovaný zdroj známé RNA/splice evidence. Samotná přítomnost záznamu
v ST2 nepřiděluje PS1 splice ani nepředvyplňuje jeho sílu. Pro pohodlnější
manuální revizi API přímo z úplné ST2 vybere řádky s multifaktoriální třídou
4 nebo 5 a zaznamenanou aberantní splice událostí. Výběr předvyplňuje jen
zdrojová fakta. Potvrzení stejné události, srovnání síly predikce a větev
Appendix J/Table 17 zůstávají povinnou manuální kontrolou.

ST2 sama nestačí k automatickému přidělení splice PS1. ENIGMA navíc vyžaduje:

- referenční P/LP klasifikaci vytvořenou podle VCEP specifications,
- přesně stejnou splice událost u hodnocené a referenční varianty,
- podobně silnou nebo silnější splice predikci u hodnocené varianty,
- správnou větev Appendix J/Table 17 podle polohy obou variant v donorovém nebo
  akceptorovém motivu,
- zahrnutí výchozího výsledku PP3 nebo PVS1 hodnocené varianty do rozhodnutí
  podle Table 17,
- u exonické varianty kontrolu souběžného proteinového následku.

Předvyplnění ze ST2 proto ponechává potvrzení stejné události, srovnání predikcí
a sílu PS1 nevyplněné.

Registr proteinových referencí:
`backend/data/ps1_protein_reference_registry.json`.
Generátor: `scripts/build_ps1_protein_reference_registry.py`.

Registr obsahuje 60 P/LP missense referencí z ST7. Aktuální sestavení obsahuje
40 záznamů `eligible` a 20 `excluded`; žádný záznam nyní není
`review_required`. Neúplný, poškozený nebo se ST7 neshodný registr zastaví start
aplikace.

Povolené klasifikační zdroje registru jsou:

- oficiální P/LP reference ENIGMA ST7 v1.2;
- verzované oficiální ENIGMA/ClinGen VCEP assertions mimo ST7;
- úplné lokální reklasifikace podle deklarované verze ENIGMA VCEP pravidel,
  jasně označené jako lokální.

Reference mimo ST7 se udržují v kurátorovaném zdrojovém souboru
`backend/data/ps1_protein_reference_extensions.json`. Generátor jej sloučí se
ST7. Validátor u oficiální assertion vyžaduje organizaci, identifikátor
assertion, verzi pravidel a datum přístupu. U lokální reklasifikace vyžaduje
posuzovatele, datum, verzi pravidel a identifikátor evidenčního záznamu.

ClinVar bez expert-panel assertion, CANVarUK, BRCA Exchange, jednotlivá studie
ani predikce samy o sobě nejsou přípustným klasifikačním základem `eligible`
reference. Mohou sloužit k nalezení kandidáta nebo jako podklady následné úplné
VCEP reklasifikace.

Identita a klasifikace současných 60 záznamů pochází ze ST7. Známá RNA evidence
se kontroluje proti úplné Table 9 a úplné ST2. SpliceAI skóre není součástí
registru. Při použití PS1 se vypočítá na požádání pro hodnocenou i referenční
variantu stejnou verzovanou službou. Transkript a normalizovaný proteinový
následek se vážou na kanonické ENIGMA RefSeq transkripty. Registr ukládá
checksum ST7, Table 9, ST2 i kurátorovaného extension souboru.

Každý záznam ukládá také podklad proteinového mechanismu. Z 40 současných
`eligible` referencí má 35 PS3 Strong funkční evidenci v Table 9. U pěti je
podkladem patogenní missense klasifikace spolu s absencí predikovaného a
potvrzeného splice efektu. Nový externí nebo lokálně reklasifikovaný záznam musí
mít odpovídající mechanismus výslovně kurátorovaný.

### 6.4 Velké exonové CNV

Runtime soubor: `backend/data/exon_cnv_evidence.json`

Zdrojový manifest: `data/sources/enigma/exon_cnv_evidence_manifest.json`

Generátor: `scripts/build_exon_cnv_evidence_snapshot.py`

Generátor stáhne celý oficiální `gnomAD-SV v2.1` BED pro GRCh37, ověří jeho
očekávanou velikost, ETag a SHA-256 a projde všechny záznamy. Seznam všech 50
exonů BRCA1/2 se generuje z úplné ENIGMA Table 4 a jejich GRCh37 kódující
intervaly z verzované souřadnicové mapy. Manifest neobsahuje c. notace variant
ani očekávaná kritéria. U každého exonu se ukládají všechny delece, které
zahrnují celý jeho kódující interval, zvlášť také shody s `FILTER=PASS`.

Za běhu se použije obecný rozhodovací graf Appendix G: typ exonová delece,
jednoznačná shoda na exon Table 4, dostupný ověřený interval, velikost nad 50 bp
a nepřítomnost jak PASS, tak filtrované zahrnující delece. Teprve při splnění
všech kroků vznikne PM2 Supporting. Přesná shoda breakpointů se nevyžaduje.

Aktuální snapshot prošel 387 477 SV záznamů, z toho 169 635 delecí. Pro BRCA2
exon 10 nebyla nalezena žádná zahrnující delece, ani filtrovaná. SHA-256
zdrojového BED je
`c843ff53b4bf36c7f733cb08565860065b3b0189375d135e33db0886381598d8`.

Soubor neobsahuje funkční kritéria. PS3/BS3 se vyhodnocují odděleně a pouze z
ENIGMA Table 9.

### 6.5 Kontrola úplnosti při startu

Table 4, Table 9, ST7, úplný ST2 splice snapshot, PS1 registry a exon-CNV
evidence jsou povinné
runtime datasety. `backend/data_validation.py` kontroluje před spuštěním API:

- existenci a čitelnost JSON,
- verzi schématu,
- očekávaný počet řádků a sloupců,
- povinná pole,
- povolené kódy a síly,
- duplicity,
- konzistenci exonových odkazů.

Neúplná nebo poškozená povinná tabulka zastaví start aplikace.

### 6.6 Referenční balík pro HGVS normalizaci

Runtime adresář je `data/reference/panel/`. Obsahuje:

- výřez oficiálního `cdot 0.2.32` RefSeq GRCh38 releasu pro přesné transkripty
  `NM_007294.4` a `NM_000059.4`, uložený beze změny cdot schématu,
- NCBI transkriptové FASTA `NM_007294.4` a `NM_000059.4`,
- nezávisle stažené NCBI proteinové FASTA `NP_009225.1` a `NP_000050.3`,
- panelový manifest, metadata, URL zdrojů a SHA-256 každého souboru.

Builder `scripts/prepare_panel_reference_bundle.py` ověří checksum celého
zdrojového cdot releasu, přesný accession s verzí, gen, protein accession,
hranice CDS, iniciační a terminační kodon a nepřítomnost interního stop kodonu.
Překlad CDS musí být totožný s proteinem uvedeným v NCBI XML i s nezávisle
staženou proteinovou FASTA. Až potom builder vytvoří malý panelový výřez a znovu
jej načte přes `cdot.JSONDataProvider`.

Celý cdot RefSeq release není načítán za běhu. Měření ukázalo přibližně 2,57 GB
RAM na worker, což není přijatelné pro čtyřworkerové nasazení. Panelový výřez má
stejné schéma a je dohledatelný k přesnému zdrojovému releasu a checksumu.

Při startu aplikace se ověří všechny checksumy, accession, transkriptová politika
a připnuté verze `hgvs` a `cdot`. Chybějící nebo změněný soubor zastaví start.
Provenance normalizace je součástí API odpovědi a ve veřejném výsledku je v
rozbalovací položce `Evidence details > Variant normalization`.

## 7. Předpočítaný coding SNV prostor

### 7.1 Soubor a rozsah

Index:

`data/precomputed/brca_module1_snv_classification_snapshot.index.json`

Metadata:

`data/precomputed/brca_module1_snv_classification_snapshot.metadata.json`

Snapshot obsahuje 47 547 coding SNV pro BRCA1 a BRCA2. Pro každou referenční pozici byly zahrnuty tři možné alternativní báze.

Rozdělení typů ve snapshotu:

| Typ | Počet |
| --- | ---: |
| initiation codon | 18 |
| missense | 35 241 |
| synonymous | 9 891 |
| nonsense | 2 397 |

Každý záznam obsahuje zejména:

- gen a `c.` notaci,
- `p.` notaci,
- typ varianty,
- GRCh37 a GRCh38 souřadnice,
- historické SpliceAI skóre pouze pro audit původního snapshotu, nikoli pro runtime klasifikaci,
- souhrn gnomAD,
- předpočítaná kritéria a třídu z okamžiku vytvoření snapshotu.

### 7.2 Použití za běhu

Runtime používá snapshot pro:

1. nezávislou kontrolu sekvenčně odvozené `p.` notace,
2. kontrolu referenční báze coding SNV,
3. lokální převod coding SNV na GRCh37 a GRCh38,
4. přístup k pomocným datům bez převzetí staré klasifikace nebo starého SpliceAI skóre.

Runtime nepřebírá předpočítanou finální třídu ani seznam kritérií jako hotový výsledek dotazu. Po kontrole vstupu se kritéria znovu vyhodnotí aktuální implementací a aktuálně načtenými runtime datasety.

Toto oddělení umožňuje použít stabilní souřadnice a nezávislou regresní kontrolu
bez zmrazení proteinového překladu nebo celé klasifikace na verzi snapshotu.

### 7.3 Stav a omezení

Metadata označují snapshot jako `snapshot_not_authoritative`. Všech 47 547
proteinových následků bylo znovu zkontrolováno připnutým lokálním HGVS enginem.
Původní zjednodušené následky 18 variant iniciačního kodonu a 15 stop-loss
variant byly nahrazeny kanonickými tvary, například `p.(Met1?)` a
`p.(Ter1864ArgextTer39)`. Reprodukovatelnou aktualizaci provádí
`scripts/refresh_snv_snapshot_protein_consequences.py`; úplný rozdílový audit
provádí `scripts/audit_hgvs_snapshot_consistency.py`.

Tento snapshot pokrývá coding SNV. Coding malé indely pokrývá samostatný snapshot popsaný níže. Exonové CNV a všechny UTR varianty pokryté nejsou.

### 7.4 Normalizovaný snapshot malých indelů

Index:

`data/precomputed/brca_normalized_indel_snapshot.index.json`

Metadata:

`data/precomputed/brca_normalized_indel_snapshot.metadata.json`

Snapshot byl vytvořen z oficiálního BRCA Exchange release 70 ze dne 8. března 2026. Tento release opravil chybná GRCh37 mapování z předchozího release. Zdrojový soubor a výsledný index mají v metadatech SHA-256.

Obsahuje 16 511 jednoznačných malých indelů. Proteinový následek `p.?` zůstává neznámý a nepoužívá se ke kontrole vstupu ani k odvození frameshift nebo in-frame typu:

| Typ | Počet |
| --- | ---: |
| frameshift | 6 898 |
| deletion, proteinový následek neznámý | 4 252 |
| insertion, proteinový následek neznámý | 1 957 |
| duplication, proteinový následek neznámý | 1 897 |
| in-frame deletion | 810 |
| in-frame delins | 269 |
| in-frame duplication | 228 |
| in-frame insertion | 113 |
| delins, proteinový následek neznámý | 87 |

Každý záznam obsahuje vstupní aliasy a kanonickou `c.` notaci, `p.` následek, typ varianty, REF a ALT alely pro GRCh37 a GRCh38, referenční transkript, BRCA Exchange release, zdrojové databáze, CA ID a VRS ID, pokud byly dostupné.

Builder `scripts/build_brca_indel_snapshot.py` čte zdroj streamovaně. Přijímá pouze BRCA1 a BRCA2 na `NM_007294.4` a `NM_000059.4`, malé indely s uvedeným proteinovým následkem, včetně `p.?`, a oběma genomovými mapováními. Konfliktní záznamy se nevkládají. Alias sdílený více záznamy se z aliasového indexu odstraní a zapíše se do metadat. Přesná kanonická notace zůstává dostupná. Aktuální release obsahuje dva takové aliasy. Runtime při startu kontroluje status, počet záznamů a checksum indexu. Chybějící nebo poškozený snapshot zastaví start aplikace.

Za běhu standardní HGVS engine ověří a kanonizuje přesný `c.` zápis přímo proti
sekvenci. Snapshot poskytuje nezávislou kontrolní `p.` notaci, aliasy používané
pro genomický reverse lookup, externí identifikátory a lokální souřadnice.
Úplný audit 16 511 indelů nalezl 8 318 záznamů se známým a shodným proteinovým
následkem a 8 193 zdrojových záznamů s neznámým následkem. Do první skupiny
patří také 15 validních variant přes hranici exonů nebo UTR, které engine
zpracuje díky HGVS normalizaci přes anotované hranice. Nebyl nalezen žádný
konflikt známých následků.
Snapshot neurčuje výslednou klinickou třídu ani automaticky nepřidává kritéria.

## 8. SpliceAI výpočet na požádání

### 8.1 Závazný profil ENIGMA Appendix J

Strojově čitelný profil je v
`data/spliceai/enigma_v1_2_spliceai_profile.json`. Hodnoty jsou převzaté z
ENIGMA BRCA1/2 VCEP v1.2 Appendix J:

- GRCh38,
- maximální vzdálenost 10 000 bází,
- unmasked režim `mask=0`,
- anotace `basic`,
- referenční transkript genu,
- maximum z `DS_AG`, `DS_AL`, `DS_DG` a `DS_DL`,
- BP4 při maximu menším nebo rovném 0,10,
- bez PP3/BP4 při hodnotě větší než 0,10 a menší než 0,20,
- PP3 při maximu větším nebo rovném 0,20.

Profil zároveň vyžaduje uložení všech čtyř delta skóre a příslušných polí
`DS_*_REF` a `DS_*_ALT`. Oficiální Appendix v původním formátu, zdrojové URL a
checksum jsou v `docs/enigma/v1.2/`.

### 8.2 Provozní zdroj a runtime cache

ARIANE nepoužívá předpočítaný prostor všech variant genu jako klasifikační
zdroj. SpliceAI se vypočítá na požádání pro konkrétní variantu přes
nakonfigurovaný ENIGMA kompatibilní endpoint. Výsledek se uloží do runtime
cache:

- `${ARIANE_RUNTIME_CACHE_DIR}/spliceai_api_cache.json`,
- na Railway `${RAILWAY_VOLUME_MOUNT_PATH}/ariane-runtime-cache/spliceai_api_cache.json`.

Stejný runtime adresář obsahuje také `bayesdel_api_cache.json` a
`coordinates_api_cache.json`. Tyto tři dynamické cache nejsou verzované v Gitu.
Načítají se před síťovým dotazem, takže dříve získaný výsledek zůstává dostupný
i při dočasném výpadku příslušného API. Klíč záznamu obsahuje identitu profilu.
Záznam s jiným profilem nebo neúplnou auditní stopou se nepoužije.

Veřejný výsledek obsahuje oddělený strukturovaný `spliceai_audit`. Hlavní klasifikace zůstává stručná a technické údaje jsou ve webovém rozhraní standardně sbalené pod položkou `Evidence details > SpliceAI`. Po rozkliknutí se zobrazí použité skóre, vybraný transkript, politika `reference_transcript`, skóre a transkript maxima přes všechny dostupné transkripty, delta pole, zdroj, GRCh38 dotaz a identifikátor cache záznamu. Stejná struktura se ukládá do auditní události dokončené klasifikace.

Zápis dynamické cache je atomický. Bez nakonfigurovaného runtime adresáře nebo
Railway volume se při lokálním vývoji používá neveřejný adresář
`.runtime-cache/` v kořeni projektu.

### 8.3 Výpočet a selhání zdroje

Cílový produkční endpoint je vlastní SpliceAI služba spuštěná z image
připnutého digestem v profilu. Veřejný Broad endpoint lze použít jako
nakonfigurovaný zdroj, není však záložním zdrojem s jinou verzí modelu. Při
změně modelu, anotace nebo referenčního genomu vznikne nový profil a nový prostor
runtime cache.

Dynamická API cache používá klíč obsahující ID závazného profilu. Staré runtime
záznamy proto nejsou znovu použity. I jednotlivá odpověď API musí výslovně
potvrdit GRCh38, vzdálenost 10 000 a `mask=0` a musí obsahovat kompletní delta,
REF a ALT pole. Jinak je skóre nedostupné. Chybějící skóre se nepřevádí na nulu
a nemůže vytvořit BP1, BP4 ani BP7.

Veřejný `spliceai_audit` je v hlavním výsledku sbalený pod
`Evidence details > SpliceAI`. Po rozkliknutí zobrazuje profil, vzdálenost,
maskování, sestavu, transkript, použitou delta hodnotu, všechna delta skóre,
REF a ALT hodnoty, zdroj, GRCh38 dotaz a cache klíč.

Interaktivní cesta používá pro Broad SpliceAI vnitřní limit 25
sekund a vnější limit 30 sekund. Pokud zdroj včas neodpoví, ARIANE dokončí
klasifikaci, označí SpliceAI jako nedostupný, nepoužije kritéria vyžadující jeho
skóre a omezení zobrazí uživateli. Požadavek tak skončí před 60sekundovým
timeoutem nginx. Obecný 12sekundový limit ostatních externích lookupů se nemění.

### 8.4 Priorita zdrojů

Používá se paměťová cache, profilově validovaná runtime cache a potom
nakonfigurovaný Broad kompatibilní výpočet. Každá
odpověď je přijata jen se stejnými parametry a úplnou auditní stopou. Neexistuje
fallback na předpočítaný prostor, starší parametry, první
dostupný transkript ani nulové skóre.

## 9. Souřadnice

### 9.1 Co znamená převod souřadnic

ARIANE potřebuje pro některé datové zdroje VCF-like alelu `chromosome`, `position`,
`REF` a `ALT` v GRCh37 nebo GRCh38. gnomAD v2.1.1 a MyVariant/BayesDel používají
GRCh37, zatímco gnomAD v3.1.2 a lokální SpliceAI data používají GRCh38.

Převod není prosté přičtení coding pozice ke genomové pozici. Musí respektovat:

- konkrétní verzi referenčního transkriptu,
- hranice exonů a intronů,
- orientaci genu; BRCA1 je na minus strand a BRCA2 na plus strand,
- normalizaci REF a ALT alely u inzercí, delecí a duplikací,
- rozdíly mezi GRCh37 a GRCh38.

Runtime proto souřadnice ručně nepočítá. Čte již normalizované souřadnice z
verzovaných dat nebo použije specializovaný HGVS resolver.

### 9.2 Lokální souřadnicové zdroje

Coding SNV mají GRCh37 a GRCh38 ve verzovaném coding SNV snapshotu. Malé indely
mají obě sestavy, REF a ALT ve verzovaném normalizovaném indel snapshotu.
Intronické SNV v podporovaném okně používají rozšířenou lokální souřadnicovou
mapu. Její referenční báze pocházejí z UCSC Genome Browser sequence API pro hg19
a hg38 a jsou navázány na stejné referenční transkripty jako ostatní snapshoty.

Při startu se intronická mapa a persistentní read-through cache načtou do
in-memory resolver cache. Aktuální implementace kontroluje tuto cache před
přímým lookupem coding SNV a indel snapshotu. Intronická mapa je při načítání
první a read-through cache nepřepisuje již existující klíč. Coding a indel
snapshot se použijí následně, pokud klíč nebyl nalezen v načtené resolver cache.

### 9.3 Síťové resolvery

Pro variantu, která nemá použitelné lokální souřadnice, je pořadí:

1. VariantValidator s referenčním transkriptem `NM_007294.4` nebo `NM_000059.4`,
2. Mutalyzer samostatně pro GRCh37 a GRCh38.

Výsledek má stav `ok`, pokud jsou dostupné obě sestavy, `partial`, pokud je
dostupná pouze jedna, a `failed`, pokud není dostupná žádná. Úspěšný síťový
výsledek se ukládá do read-through cache. Úplné přechodné selhání se necachuje,
aby mohl pozdější dotaz resolver zopakovat.

Nejednoznačný výsledek se nesmí převést na první nalezené ID. Exonová CNV s
neurčitými breakpointy síťové resolvery vůbec nevolá, protože nemá jednu
konkrétní VCF alelu.

Souřadnice jsou verzovaná vlastnost kombinace referenčního transkriptu,
genomového sestavení a normalizační politiky. Nemění se při každém dotazu, ale
musí se znovu validovat při změně kterékoliv z těchto částí.

## 10. gnomAD data

Runtime soubory jsou v `backend/data/gnomad/`.

- `gnomad_brca_frequency_snapshot.json`: frekvenční záznamy pro panel,
- `gnomad_brca_coverage_snapshot.json`: per-position coverage,
- `gnomad_panel_manifest.json`: verze zdrojů, intervaly a klasifikační politika.

Používané datasety a frekvenční data:

- gnomAD v2.1.1 exomes non-cancer na GRCh37,
- gnomAD v3.1.2 genomes non-cancer na GRCh38,
- samostatný per-position coverage snapshot.

Aktivní zdroje, jejich verze, GCS identity, panelové intervaly a genově
specifické frekvenční politiky jsou v
`backend/data/gnomad/gnomad_panel_manifest.json`. Intervaly ani prahy BA1, BS1
a PM2 nejsou zapsané jako globální konstanty v runtime. Samotné přidání
intervalu nesmí aktivovat nový gen.

Snapshot se sestavuje skriptem
`scripts/refresh_gnomad_panel_snapshot.py` přímo z oficiálních gnomAD Hail
Tables. Hail načte jen intervaly panelu včetně definovaného paddingu. Celý
genomový balík se nemusí kopírovat na lokální disk. Build používá anonymní
read-only přístup k veřejnému GCS bucketu a připnutou verzi Hail z
`requirements-data.txt`.

gnomAD v2.1.1 poskytuje non-cancer FAF95 přímo v poli `faf` oficiální Hail
Table. Builder načte hodnoty pro `afr`, `amr`, `eas`, `nfe` a `sas` a uloží
jejich maximum.

gnomAD v3.1.2 poskytuje pro non-cancer subset populační AC a AN, ale neposkytuje
hotové subsetové FAF95. Builder proto používá přímo
`hail.experimental.filtering_allele_frequency` nad oficiálními non-cancer AC a
AN. Výpočet probíhá uvnitř Hail, nikoli vlastní numerickou implementací v
ARIANE. Použije se jednostranný 95% Poissonův interval a pět ENIGMA
non-founder populací.

FAF95 vrácené gnomAD Browser API není subsetově specifické. API vrací stejnou
hodnotu pro hlavní, non-cancer i controls-and-biobanks dataset. Tato hodnota se
proto nepoužívá jako non-cancer FAF95. API se nepoužívá ani jako frekvenční
fallback.

Coverage se při stejném buildu obnovuje z oficiálních Hail coverage Tables pro
obě sestavy. Pro v2 se používá exome coverage release 2.1. Pro v3 se používá
genome coverage release 3.0.1, který používá také oficiální datová pipeline
gnomAD Browseru pro dataset gnomAD v3.

Každý záznam obsahuje FAF95 po populacích, maximum, populaci maxima, scope a
metodu. Obsahuje také AC a AN pěti skórovaných non-founder populací, příznak
jejich skutečné přítomnosti a oddělený kontext vyloučených skupin. V aktuálním
panelovém výřezu je 1 644 záznamů pozorovaných pouze ve vyloučených skupinách.
Metadata obsahují Hail URI, URL metadata objektu, ETag, GCS MD5/CRC32C,
generation, velikost, datum zdroje, počet záznamů a zákaz fallbacku na raw AF.
Frekvenční i coverage JSON obsahují SHA-256 kanonického obsahu a SHA-256
manifestu. Runtime při načtení kontroluje oba datasety, povolenou metodu FAF95,
scope, identity zdrojů, manifest a checksumy. Při selhání se cache nepoužije a
degradace se zobrazí uživateli.

Kontrola dostupných release nic sama neaktivuje:

```bash
python scripts/refresh_gnomad_panel_snapshot.py check-updates
```

Nový release vyžaduje kontrolu souladu s pravidly ENIGMA, porovnání výsledků a
schválení regresních testů. Aktivní release se nemění pouze proto, že se v
bucketu objevil novější adresář.
Kontrola navíc sestaví cestu ke stejnému small-variant sites Hail Table produktu
a ověří jeho dostupnost. gnomAD 3.1.3 je například release tandemových repetic,
nikoli náhrada small-variant dat 3.1.2, a proto se jako použitelná aktualizace
nehlásí.

Obnova v samostatném datovém prostředí:

```bash
pip install -r requirements-data.txt
python scripts/refresh_gnomad_panel_snapshot.py refresh
python scripts/refresh_gnomad_panel_snapshot.py validate
```

Obnova pomocí připnutého Docker image:

```bash
docker run --rm -v "$PWD:/workspace" -w /workspace \
  hailgenetics/hail:0.2.137 \
  python3 scripts/refresh_gnomad_panel_snapshot.py refresh
```

Builder nejprve vytvoří pracovní exporty. Runtime JSON nahradí atomicky až po
úspěšné validaci obou frekvenčních datasetů, obou coverage datasetů a checksumů.
Starý `scripts/build_gnomad_v3_brca_snapshot.py` je ponechán jen pro audit
předchozí numerické implementace a nesmí publikovat runtime data.

BA1 a BS1 vyžadují průměrnou hloubku alespoň 20. PM2 vyžaduje prokázanou
nepřítomnost v obou požadovaných non-cancer datasetech a průměrnou hloubku
alespoň 25. Samotná absence varianty v JSON není dostačující. PM2 se podle
ENIGMA v1.2 nepoužívá pro indely.

Coverage lookup agreguje všechny pozice genomového rozsahu `REF`. ENIGMA v1.2
neuvádí šířku flanking okna, proto se okolní báze mimo vlastní variantu
nepřidávají. Tato definice je uvedena v diagnostických polích jako
`coverage_scope: variant_reference_span`. Pro SNV je výsledek shodný s mean
depth na lokusu, který zobrazuje gnomAD. Pro vícebázový rozsah musí být v cache
všechny pozice, jinak se vyhodnocení uzavře bez PM2.

### 10.1 Otevřené implementační body

1. ENIGMA v1.2 pro PM2 požaduje průměrnou hloubku alespoň 25 v oblasti kolem
   varianty, ale neurčuje šířku této oblasti. Současná reprodukovatelná definice
   ARIANE používá genomový rozsah alely `REF`. Pro SNV je to jedna pozice.
   Případná změna na flanking okno se nesmí provést odhadem. Vyžaduje potvrzení
   ENIGMA, verzovanou změnu politiky a nové regresní testy.

2. ENIGMA zakazuje BA1 a BS1 u dobře doložených patogenních founder variant,
   ale neposkytuje jejich úplný strojově čitelný seznam. Lokální kurátorovaný
   soubor proto obsahuje pouze varianty s doloženým zdrojem. Nepřítomnost
   varianty v tomto souboru sama o sobě nedokazuje, že nejde o founder variantu.
   Rozšíření seznamu musí obsahovat kanonickou HGVS notaci, používaný transkript,
   founder populaci, tvrzení o patogenitě, zdroj, datum přístupu a checksum.

### 10.2 Kompatibilita coverage pro gnomAD v3

Frekvenční data ARIANE pocházejí z gnomAD v3.1.2. Oficiálně veřejně dostupná
genomová coverage tabulka pro GRCh38 je však označena jako gnomAD r3.0.1.
gnomAD v3.1 následně přidal 4 598 vstupních gVCF a výsledný callset obsahuje
76 156 genomů.
Release v3.1.2 opravil některé genotypy a frekvence, ale samostatnou coverage
tabulku v3.1.2 neuvolnil.

Použitá r3.0.1 coverage je oficiální gnomAD produkt a jediná veřejně
distribuovaná genome coverage tabulka pro řadu v3. Není však vzorkově totožná s
variantním callsetem v3.1.2. ENIGMA v1.2 požaduje gnomAD v3.1 a pokrytí oblasti
kolem varianty, ale neuvádí konkrétní coverage release. Z toho nelze odvodit,
že ENIGMA rozdíl výslovně schválila.

Do získání potvrzení od ENIGMA je stav veden jako doložený provozní zdroj s
nepotvrzenou přesnou kompatibilitou. Zdroj a verze se vždy zobrazují v auditu a
nesmí být přepsány na `v3.1.2 coverage`. Podklady:

- [gnomAD v3.1 release](https://gnomad.broadinstitute.org/news/2020-10-gnomad-v3-1-new-content-methods-annotations-and-data-availability/),
- [gnomAD v3.1.2 minor release](https://gnomad.broadinstitute.org/news/2021-10-gnomad-v3-1-2-minor-release/),
- `gs://gcp-public-data--gnomad/release/3.0.1/coverage/genomes/gnomad.genomes.r3.0.1.coverage.ht`.

### 10.3 Přidání dalších genů

Builder gnomAD umí načíst libovolný počet genových intervalů z položky
`targets` v manifestu. To samo o sobě nestačí k povolení klasifikace dalšího
genu. Manifest používá schéma 2. Každá položka `targets` musí explicitně
odkazovat přes `classification_policy_id` na aktivní politiku. Žádná politika
se nedědí automaticky a samotný interval nestačí. BRCA1 a BRCA2 nyní explicitně
odkazují na `enigma_brca_v1_2`. Runtime čte z této politiky prahy, síly, body,
požadované datasety, populační skupiny, pokrytí a vyloučené typy variant.
Rozsah této politiky je výslovně `gnomad_frequency_criteria_only`: manifest
definuje BA1, BS1 a PM2, nikoli ostatní VCEP kritéria. Ta musí mít pro nový gen
samostatnou implementaci a validaci a z BRCA1/2 se nikdy nedědí.
Politika navíc obsahuje uzavřený seznam `applicable_genes`. Pokus připojit
například nový gen k `enigma_brca_v1_2` builder i runtime odmítnou.

Politika dalšího genu musí určit:

- VCEP a přesnou verzi pravidel,
- referenční transkript s verzí,
- intervaly GRCh37 a GRCh38,
- povolené gnomAD release, callsety a populační skupiny,
- prahy BA1, BS1 a PM2 včetně požadavků na pokrytí,
- pravidla pro founder varianty a zdroj jejich kurátorovaného seznamu,
- kritéria, která se pro daný gen nesmějí automatizovat.

Builder odmítne cíl bez aktivní politiky, bez referenčního transkriptu nebo bez
provenance VCEP. Runtime pro neznámý gen nebo neshodnou politiku vrátí stav
`policy_unavailable` a nepoužije BA1, BS1 ani PM2. BRCA politika se na nový gen
nikdy nepřenese jako fallback.

Nový gen se nesmí aktivovat pouhým přidáním intervalu. Aktivace vyžaduje nový
datový build, validační matici pro daný gen, kontrolu všech gene-specific
tabulek a explicitní povolení genu v API a uživatelském rozhraní. Do té doby
zůstávají jedinými povolenými geny BRCA1 a BRCA2.

### 10.4 Výsledek auditu gnomAD proti ENIGMA v1.2

| Oblast | Stav |
| --- | --- |
| gnomAD v2.1.1 exomes non-cancer | Oficiální Hail data, připnutý zdroj a identita objektu |
| gnomAD v3.1.2 genomes non-cancer | Oficiální Hail data, připnutý zdroj a identita objektu |
| FAF95 pro v2 | Oficiální hodnoty z Hail tabulky |
| FAF95 pro v3 non-cancer | Výpočet oficiální Hail funkcí z oficiálních AC a AN, bez náhradní metriky |
| Populace pro skórování | AFR, AMR, EAS, NFE a SAS podle Appendix G |
| Founder a ostatní vyloučené populace | Uloženy odděleně jako auditní kontext, nevstupují do FAF95 maxima |
| BA1 a BS1 | Prahy, QC, hloubka a founder výjimka jsou v genově specifické politice |
| PM2 | Oba požadované datasety, hloubka 25, pouze povolené typy variant |
| Jediné outbred pozorování | Neinformativní, nevytváří BA1, BS1 ani PM2 |
| Raw AF a popmax AF | Pouze auditní kontext, nikdy fallback pro BA1 nebo BS1 |
| Runtime síťový fallback | Není implementován |
| Neúplný nebo starý snapshot | Odmítnut podle manifestu, provenance a checksumu |

Audit neumožňuje tvrdit, že jsou vyřešeny všechny interpretační nejasnosti.
Zůstávají tři explicitní omezení uvedená také přímo v manifestu: ENIGMA
neurčuje šířku coverage okna, veřejná v3 coverage je z release 3.0.1 a seznam
patogenních founder variant není v ENIGMA publikován jako úplný strojově
čitelný katalog. ARIANE tyto mezery nenahrazuje odhadem ani skrytým fallbackem.

## 11. BayesDel a AlphaMissense

BayesDel_noAF a AlphaMissense se získávají jedním dotazem nad genomovou variantou přes MyVariant.info a ukládají se do lokální cache.

Do trvalé cache se ukládá pouze úspěšná anotace nebo explicitní odpověď, že
varianta v MyVariant nebyla nalezena. Odpověď `no_score`, při které služba
variantu vrátí bez BayesDel i AlphaMissense anotace, se považuje za opakovatelný
neúplný výsledek. Neukládá se a následující požadavek provede nový dotaz. Při
načtení cache se ignorují také starší prázdné a `no_score` záznamy.

BayesDel se používá pouze v rozhodovacích větvích PP3 a BP4 popsaných výše. AlphaMissense se vrací jako doplňující anotace a samo o sobě nevytváří samostatné ENIGMA kritérium.

Selhání služby, chybějící GRCh37 souřadnice nebo nenalezená anotace mají odlišné stavové kódy. Důvod se zobrazí v diagnostice. Chybějící BayesDel nemůže být nahrazen předpokládanou hodnotou.

## 12. ClinVar a ClinGen

ClinVar a ClinGen ERepo se používají pro externí srovnání, auditní kontext a
předvyplnění ověřitelných faktů v manuální revizi proteinového PS1. Jejich
klasifikace se automaticky nepřičítá jako ACMG nebo ENIGMA kritérium.

Pokud vyhledávání vrátí více kandidátů bez jednoznačné přesné shody, stav je `ambiguous`. Aplikace nevybere první ID.

## 13. Degradované datové zdroje a fail-closed chování

Povinné ENIGMA tabulky zastaví start aplikace, pokud jsou neúplné nebo nečitelné.

U ostatních zdrojů se degradace ukládá do centrálního registru a zobrazuje v
klasifikačním výsledku a `/api/health`. Hlášení centrálního registru obsahuje
komponentu a konkrétní důvod. Absolutní serverové cesty se pro uživatele
zkracují na tvar `…ariane/...`.

Podrobné odpovědi síťových resolverů, HTTP těla a výjimky se zapisují do
serverového logu. Veřejný klasifikační výsledek je neopakuje. Uživatel dostane
jedno stručné a akční shrnutí pro nedostupné souřadnice, ClinVar nebo ClinGen.
Tím se zachová auditovatelná technická příčina bez zveřejnění interních cest a
dlouhých odpovědí poskytovatelů.

Bez potvrzeného vstupu nebo bez dat požadovaných konkrétním pravidlem se dané kritérium nepoužije. Aplikace nesmí nahrazovat chybějící data nulou, fixture hodnotou, prvním nalezeným ID nebo ručním fallback slovníkem.

## 14. Audit a reprodukovatelnost

Každý předpočítaný dataset má být doprovázen metadaty obsahujícími alespoň:

- název a verzi datasetu,
- datum vytvoření,
- zdroj,
- referenční transkripty a genomová sestavení,
- počet záznamů,
- checksum,
- stav úplnosti,
- verzi nebo commit generátoru, pokud je k dispozici.

Změna oficiální ENIGMA verze, transkriptu, genomového sestavení, predikčního modelu nebo klasifikační logiky vyžaduje nové sestavení dotčených snapshotů a opakování regresních testů.

### 14.1 Regresní matice typů variant

Soubor `tests/test_variant_type_regression_matrix.py` udržuje explicitní regresní matici podporovaných typů variant. První část ověřuje odvození normalizovaného typu pro nonsense, frameshift, missense, synonymous, intronickou SNV, canonical splice-site variantu, in-frame indel, exonovou deleci a duplikaci, iniciační kodon a obě UTR oblasti.

Druhá část pro reprezentativní varianty kontroluje celý deterministický průchod Module 1 od typu varianty přes lokální tabulky a snapshoty po aplikovaná kritéria a výslednou třídu. Každý řádek obsahuje očekávaná i zakázaná kritéria. Tím se kontroluje nejen přítomnost správné větve, ale také nepřípustné přičtení kritéria z jiné větve, například PP3 u nonsense nebo frameshift varianty.

Synonymous varianty mají v matici obě ENIGMA větve. Varianta uvnitř funkční domény s potvrzeným SpliceAI nejvýše 0,1 vyžaduje BP4 a BP7 a nesmí dostat BP1. Varianta mimo funkční domény se stejným limitem SpliceAI vyžaduje BP1 Strong a nesmí současně dostat BP4 ani BP7. Pozitivním regresním případem pro BP1 je `BRCA1 c.306A>C p.(Ala102=)`.

RNA-dependent canonical splice-site a initiation-codon větve samostatně ověřují, že aplikace nevytvoří automatické PVS1 a místo něj vrátí doporučení ke strukturované manuální revizi. Matice nepoužívá živé externí služby. Používá verzované lokální ENIGMA tabulky a snapshoty, takže musí být deterministická.

### 14.2 Interakce funkční, RNA a predikční evidence

ARIANE nepoužívá plošné potlačení podle příznaku `has_functional_evidence`. Interakce se vyhodnocují podle mechanismu a hierarchie důkazů z ENIGMA v1.2 Figure 1A, Figure 1B, Figure 1C a Appendix E.

Přijaté PVS1 (RNA) nahrazuje slabší bioinformatické kódy pro stejný experimentálně potvrzený splice důsledek. Přijaté BP7 Strong (RNA) nahrazuje BP7 Supporting, ale podle Figure 1B obecně zachovává ostatní použitelné bioinformatické kódy. Před přijetím BP7 Strong (RNA) se samostatně kontroluje typ varianty a funkční doména. Missense varianta uvnitř ENIGMA funkční domény musí mít v původním automatickém výsledku aplikované BS3 s Table 9 provenance. Ruční potvrzení tuto podmínku nemůže nahradit. PS3 nebo BS3 bez PVS1 automaticky nepotlačuje PP3, BP4, BP7 ani BP1, protože Figure 1C výslovně požaduje zachování relevantních bioinformatických kódů.

Každé nahrazení, zachovaná potenciální interakce nebo konflikt se vrací ve strukturovaném poli `evidence_interactions`. Ve webovém rozhraní je zobrazeno v rozbalovací části `Evidence interaction warnings`. Přesná matice je v `docs/evidence_interaction_matrix.md`.

Appendix Table 11 se používá jen pro variantově specifické klinické anotace.
Aktuálně z ní ARIANE deterministicky vybere záznam, který ENIGMA výslovně
označuje jako `Proven reduced penetrance allele`, tedy BRCA1 c.5096G>A
p.(Arg1699Gln). Ve výsledku se zobrazí `Variant with reduced penetrance`, odkaz
na Appendix Table 11 a publikace PMID 22889855 a 28490613. Anotace není
kritérium, nepřidává body a nemění třídu. Pokud by zdrojová tabulka změnila
sloupce nebo explicitní označení zmizelo, startup validace selže.

Související dokumenty:

- `docs/enigma_source_data_audit.md`,
- `docs/evidence_interaction_matrix.md`,
- `docs/manual_evidence_review.md`,
- `docs/vus_explanation_and_golden_cases.md`.

## 15. Architektura spuštění klasifikace

ARIANE používá typovaný DAG, který při chybě nevytvoří klasifikaci. Normalizace
vstupu probíhá samostatně před vyhodnocením evidence. Získání dat, rozhodnutí o
kritériích, řešení interakcí evidence, výsledná klasifikace a prezentace jsou
oddělené vrstvy.

Produkční klasifikaci provádí provider graf `ariane.vcep.classification`, verze
`4.0.0-gene-policy-provider-dag`. Starý sekvenční evaluator se v aplikační cestě neimportuje
ani nespouští. Jediná povolená hodnota `ARIANE_CLASSIFIER_ENGINE` je `dag` a jde
zároveň o výchozí hodnotu. Režimy `legacy`, `shadow` ani fallback nejsou
dostupné.

DAG používá typy `NormalizedVariant`, `EvidenceBundle`,
`CriterionFamilyResult`, `CriterionDecision` a `VariantAssertion`. Jednotlivé
rodiny pravidel nevkládají výsledky do společně mutovaného slovníku. Každá
vrací samostatný neměnný výsledek. Sloučení probíhá v uzlu pro interakce
evidence.

### 15.1 Produkční automatický graf

| Uzel | Úloha |
| --- | --- |
| `contract.classification_request` | Kontrola normalizovaného požadavku |
| `provider.coordinates` | Souřadnice GRCh37 a GRCh38 |
| `provider.spliceai` | SpliceAI pro schválený referenční transkript |
| `provider.bayesdel` | BayesDel_noAF a informační AlphaMissense |
| `provider.gnomad` | Populační frekvence a pokrytí |
| `provider.enigma.table9` | Lookup ve validovaném Table 9 datasetu, záznam verze a checksumu |
| `provider.clinical_lr`, `provider.exon_cnv`, `provider.protein_ps1` | Další klasifikační evidence |
| `contract.evidence_bundle` | Sestavení a kontrola typované evidence |
| `context.spliceai.provenance` | Použití konfigurovaného SpliceAI skóre a auditní porovnání s Table 9 bez přepsání výsledku |
| `rule.population_frequency` | BA1, BS1 a PM2 z předané gnomAD evidence |
| `rule.exon_cnv.population` | Populační větev pro exonové delece a duplikace |
| `rule.functional.table9` | PS3 nebo BS3 a Figure 1C decision path |
| `rule.pvs1_pm5` | PVS1, PVS1 RNA a PM5 PTC |
| `rule.clinical_lr` | PP4 nebo BP5 z validované kombinované klinické LR evidence |
| `rule.protein_ps1` | Proteinové PS1 ze schváleného registru referencí |
| `rule.bioinformatic.figure1a` | PP3, BP4, BP7 a BP1 podle Figure 1A |
| `policy.evidence_interactions` | Mechanism-aware zachování nebo potlačení překrývající se evidence |
| `policy.enigma_combination` | BA1 stand-alone, ENIGMA Table 3 nebo bodová metoda pro mixed evidence |
| `review.manual_triage` | RNA, splice PS1, protein PS1 a initiation review výstupy |
| `contract.variant_assertion` | Kontrola úplnosti a vnitřní konzistence výsledku |
| `projection.public_result` | Převod typovaného výsledku do veřejného API kontraktu |

### 15.2 Graf pro ručně doplněnou evidenci

Přepočet po ručně doplněné odborné evidenci používá samostatný graf
`ariane.vcep.manual-evidence`, verze `2.0.0-gene-policy-dag`.

| Uzel | Úloha |
| --- | --- |
| `contract.manual_evidence_inputs` | Kontrola struktury automatických a ručních kritérií |
| `rule.manual_evidence` | Kontrola povinných podkladů a výpočet povolené síly |
| `policy.manual_evidence_interactions` | RNA a proteinová deduplikace podle mechanismu |
| `policy.manual_enigma_combination` | Review-adjusted ENIGMA klasifikace |

Frontend ruční evidence neobsahuje kopii prahů ani funkci pro odvození síly.
Odesílá původní automatická kritéria, kontext varianty a surová pole formuláře
do backendového grafu. Způsobilost kritéria, úplnost podkladů, sílu, body a
interakce evidence vrací výhradně backend. Také úplnost formuláře, tedy výběr
kritéria, jméno hodnotitele, datum, poznámku a zdroj, kontroluje backendový API
model. Klient pouze odesílá surová pole a zobrazuje odpověď.
Auditní export přebírá odvozené síly z backendového `amended_working_result` a
nevytváří vlastní výpočet.

Pokud automatický výsledek obsahuje BA1, zůstává v manuálním grafu
stand-alone benigní klasifikací. Ručně doplněná kritéria se zachovají v
auditní stopě, ale BA1 nepřepínají do bodové klasifikace mixed evidence.

Neúplný odborný podklad vrací validační chybu 422. Neočekávaná interní chyba
uzlu vrací 503. V obou případech se nová klasifikace nevydá a audit obsahuje
identifikátor chybného uzlu. Chybějící podklad se nenahrazuje odhadem, nulou,
fixture hodnotou ani fallbackem.

### 15.3 Hranice získávání dat

Souřadnice, Table 9, gnomAD, SpliceAI, BayesDel, klinické LR, PS1, důležitá
rezidua a exonové CNV získávají samostatné provider uzly uvnitř produkčního
DAGu. Každý provider vrací hodnotu, stav dostupnosti a provenance. Teprve uzel
`contract.evidence_bundle` sestaví typovanou evidenci pro pravidlové uzly.

Síťové API, lokální cache a předpočítané datasety zůstávají zdroji providerů,
nikoli součástí pravidlových uzlů. Nedostupnost se nesmí převést na nulové
skóre, nesplněné kritérium, první nalezené ID nebo náhradní fixture.

### 15.4 Starý klasifikátor

Produkční kód `backend/main.py` neimportuje `backend/modules/classifier.py`.
Starý klasifikátor zůstává pouze jako nezávislý testovací oracle pro kontrolu
parity. Jeho fyzické odstranění vyžaduje přesun zbývajících testů na veřejný
DAG kontrakt a samostatné schválení regresní parity. Odstranění nesmí být
spojeno se změnou klinických pravidel nebo očekávaných výsledků.

### 15.5 Regresní ověření DAG

Paritní testy pokrývají tutorialové varianty, hlavní typy variant, PVS1 a PM5,
Table 9, proteinové PS1, exonové CNV, klinické LR, BA1 terminální větev, BS1
mixed evidence, PM2 a manuální RNA interakce. Aktuální počet testů je uváděn v
protokolu konkrétního vydání, ne jako vlastnost architektury.

Podrobný návrh, invarianty a plán odstranění starého souboru jsou v
`docs/classification_dag_architecture.md`.

### 15.6 Veřejný přístup k pravidlům a tabulkám

ARIANE používá verzovaný katalog `backend/data/enigma_rule_catalog.json` jako
společný registr pravidel, oficiálních zdrojů a jejich checksumů. Souřadnice
uzlů a spojnic dalších diagramů jsou odděleny v
`backend/data/enigma_rule_diagrams.json`. Prezentační geometrie je tím oddělena
od klasifikační logiky a lze ji upravit bez změny výpočtu kritérií.

Diagramy mají povinnou provenienci:

- `official_redraw` označuje přístupné SVG překreslení konkrétní ENIGMA Figure;
  záznam vždy odkazuje na původní obrazový panel a oficiální dokument;
- `ariane_derived` označuje rozhodovací cestu sestavenou ARIANE z uvedené
  tabulky, appendixu nebo části specifikace. Nejde o novou ENIGMA Figure a
  rozhraní ji tak nesmí označovat.

Překreslené oficiální podklady pokrývají Figure 1A, 1B a 1C ze Specifications,
Appendix Figures 3 až 6 pro PVS1, exonové a PTC/PM5 mapy z Appendix Figures 1,
2, 7 a 8 a mechanism-aware strom Appendix Figure 9. Diagramy bez samostatné
oficiální předlohy pokrývají populační evidenci, způsobilost proteinové PS1
reference, combined clinical LR pro PP4/BP5 a finální volbu klasifikačního
postupu. Každý je v rozhraní viditelně označen jako odvozený diagram ARIANE.

Specializovaný endpoint Table 9 čte stejný validovaný runtime dataset jako
klasifikátor. Obecný prohlížeč navíc zpřístupňuje všech 42 oficiálních tabulek z
odděleného prezentačního snapshotu.

Stránka `ENIGMA rules` zobrazuje rozhodovací stromy, všechny Specification,
Appendix a Supplementary Tables, původní obrázky a verze zdrojů. U aplikovaných
PVS1 pro PTC, kanonické splice varianty, exonové
delece a duplikace, PVS1 (RNA), PS3, BS3, PP3, BP4, BP7 a BP1 lze
rozbalit skutečnou rozhodovací cestu vytvořenou během výpočtu. Ostatní diagramy
zatím slouží jako referenční zobrazení bez předstírané runtime stopy. API vrací
pouze veřejná pole, stránkuje každou tabulku a nikdy nezveřejňuje
lokální cesty serveru. Endpointy a datový
kontrakt jsou popsány v `docs/classification_dag_architecture.md`.

Rozbalená rozhodovací cesta je v tabulce aplikovaných kritérií umístěna do
samostatného panelu přes plnou šířku tabulky. Lineární cesta používá nejvýše tři
uzly na řádek a pokračuje střídavě na dalším řádku. Panel proto nepotřebuje
vodorovný posuvník a zachovává celý text uzlů i pozorovaných hodnot.

Rozhodovací diagramy se zobrazují jako responzivní SVG bez vodorovného posuvníku.
Oficiální Figure používají ručně zadanou geometrii podle předlohy. Každý uzel
může mít vlastní rozměr a každá spojnice vlastní body a pozici popisku. Použitá
cesta obsahuje hodnoty, které podpořily jednotlivé volby, a je barevně odlišena
od nepoužitých větví, pokud klasifikační modul vrací strukturovaný
`decision_path`. Diagram bez připojené runtime stopy se zobrazuje jako obecná
referenční cesta a nesmí předstírat, že byla pro aktuální variantu použita.

SVG renderer vede spojnice jako tenké pravoúhlé trasy se zaoblenými rohy.
Popisky větví mají vlastní malé pozadí a skutečné hodnoty rozhodnutí jsou v
samostatných víceřádkových anotacích nad uzlem, aby se nepřekrývaly s otázkou v
uzlu. Anotace nesmí text zkracovat. Výška anotace a horní prostor SVG se počítají
z celého textu runtime rozhodnutí.

Samostatná galerie obsahuje původní obrazové panely z checksumem ověřených kopií
ENIGMA Specifications a Appendix v1.2. Appendixová média se reprodukovatelně
extrahují skriptem `scripts/extract_enigma_appendix_figures.py`; skript pouze
přejmenuje embedded média a nijak nemění jejich obsah.
