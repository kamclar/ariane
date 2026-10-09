# Dohled nad aktualitou zdrojů

ARIANE odděluje tři různé otázky, které se nesmějí zaměňovat:

1. **Integrita aktivních dat:** startup a `/api/health` ověřují, že lokální
   snapshoty existují, odpovídají manifestům, checksumům a očekávanému schématu.
2. **Aktualita upstream zdroje:** read-only source watch zjišťuje, zda oficiální
   zdroj zveřejnil novou verzi nebo změnil již sledovaný záznam.
3. **Aktivace nové verze:** samostatný, odborně schválený rebuild rozhodne, zda je
   kandidát kompatibilní s aktivní ENIGMA politikou a může vstoupit do produkce.

Nedostupnost externího webu proto nezastaví běžnou klasifikaci. Runtime vždy
používá poslední schválený a checksumem ověřený lokální snapshot. Source watch
nemá právo měnit produkční soubory ani manifesty.

## Centrální registr

Soubor `data/sources/source_watch_registry.json` uvádí sledované zdroje,
oficiální endpoint, frekvenci kontroly, aktivní baseline a proces pro vytvoření
kandidáta. Registr výslovně vyžaduje:

- `automatic_activation: false`,
- `runtime_network_dependency: false`,
- archivaci identity kandidátního zdroje,
- datový a kriteriální diff,
- cílené i úplné regresní testy,
- zaznamenané odborné schválení před aktivací.

Automaticky lze rozumně kontrolovat:

| Zdroj | Co se porovnává | Frekvence |
|---|---|---|
| ENIGMA BRCA1/2 VCEP CSpec | verze, released stav, revize a datum změny GN092/GN097 | týdně |
| ClinGen ERepo ENIGMA VCEP | fingerprint ENIGMA BRCA1/2 řádků, audit SHA-256 úplného exportu a verze Summary API | týdně |
| UCSC ENIGMA clinical LR | SHA-256 aktuálního BRCAmfa BigBed a `dataVersion` | měsíčně |
| gnomAD | identita aktivních Hail produktů a dostupnost stejného produktu v novějších release | měsíčně |

SpliceAI profil, referenční transkripty, souřadnicové zdroje, founder varianty,
klinicky významná rezidua a exon-CNV registr jsou ve stejném reportu vedeny jako
`manual_only`. Nemají jeden bezpečný strojově čitelný endpoint, jehož nejvyšší
číslo verze by bylo možné bez metodického posouzení označit za vhodnou
aktualizaci. Jejich kontrolu spouští změna VCEP pravidel, referenčního
transkriptu, genomové sestavy, modelu, licence nebo nová kurátorská publikace.

## Spuštění a interpretace

Kompletní read-only kontrola:

```bash
python scripts/check_source_updates.py \
  --output source-watch-report.json \
  --markdown-output source-watch-report.md
```

Pouze jeden zdroj, například gnomAD:

```bash
python scripts/check_source_updates.py --source gnomad_population_data
```

Lokální kontrola registru bez sítě:

```bash
python scripts/check_source_updates.py --validate-registry
```

Návratové kódy:

| Kód | Význam |
|---:|---|
| `0` | kontrolované automatické zdroje odpovídají baseline; manuální zdroje jsou pouze evidovány |
| `2` | byla nalezena nová verze, změněná upstream identita nebo je splatná manuální revize |
| `3` | alespoň jeden upstream zdroj nebylo možné zkontrolovat |
| `4` | chyba nebo nebezpečné nastavení centrálního registru |

Stav `check_failed` neznamená, že aktivní data jsou chybná. Znamená jen, že
aktuálnost nebylo možné v daný okamžik potvrdit. Kontrola se zopakuje; aktivní
snapshot zůstává beze změny.

GitHub workflow `.github/workflows/source-watch.yml` spouští kontrolu každé
pondělí a lze jej spustit ručně před releasem. Výsledek zapisuje do job summary.
Nález změny i technicky neúspěšná kontrola workflow zviditelní jako neúspěšný
job, ale nic neaktivuje a nepotřebuje přístup k produkčnímu serveru.

## Postup po nálezu změny

1. Ověřit, že nejde jen o technickou změnu endpointu nebo transientní chybu.
2. Uložit kandidátní zdroj, čas stažení, URL, verzi a checksum mimo aktivní data.
3. Přečíst oficiální release notes a posoudit vztah k ENIGMA VCEP v1.2.
4. Sestavit kandidátní snapshot a vytvořit diff přidaných, odebraných a změněných
   záznamů, prahových přechodů a změn přidělených kritérií.
5. Provést odbornou revizi, cílené regrese a kompletní test suite.
6. Teprve poté v jednom reviewovatelném commitu změnit zdroj, manifest, odvozená
   data, testy, verzi aplikace a dokumentaci.

Novější neznamená automaticky vhodnější. To je důležité zejména u gnomAD: jiný
major release může být podle ENIGMA přípustný, ale má jiné kohorty, coverage a
datové schéma. Hlídač jej uvede jako informaci; migrace na gnomAD v4 je
samostatná validační změna, nikoli automatický update v2/v3 snapshotů.
