## Beispiele für Kapitel 4 - Petrov-Galerkin-Verfahren

In den Unterordnern `example{1,2,3}` befinden sich die drei am Ende von Kapitel 4 betrachteten Beispiele für die Petrov-Galerkin-Diskretisierung. 
Weiter finden sich in `shared` mehrere Skripte, die für alle drei Beispiele benötigt werden.

Die Beispiele enthalten neben den Startskripten `ch4ex{1,2,3}_run.m` jeweils auch die bereits vorberechneten Daten `ch4ex{1,2,3}_precomputed.mat`, da die benötigten Berechnungen 30 Minuten und mehr dauern können.

Um alle Beispiele nacheinander auszuführen, empfiehlt es sich, das Skript `shared/ch4_precompute_examples.m` zu starten.
Dieses führt die benötigten Berechnungen aus und speichert nur für die Auswertung relevanten Daten (bei der Berechnung werden knapp 3-4 GB Arbeitsspeicher benötigt, da die Eigenwertprobleme sehr groß werden).

Die in Kapitel 4 vorkommenden Abbildungen können mit dem Skript `shared/ch4_load_precomputed_data.m` nachgestellt werden. Um weiter die zugehörigen `tikz`-Dateien zu generieren, muss die Variable `generate_tikz` auf den Wert `true` gesetzt werden.
