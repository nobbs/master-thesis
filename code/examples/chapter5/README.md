## Beispiele für Kapitel 5 - Reduzierte-Basis-Methode

Die Unterordner `example{1,2}` enthalten die beiden Beispiele, die am Ende von Kapitel 5 zu finden sind. Neben den benötigten Skripten um die Berechnungen auszuführen finden sich auch bereits vorberechnete Daten, da die Berechnungen vor allem bei zweitem Beispiele mehrere Stunden dauern.
Beide Beispiele enthalten jeweils eine Datei vom Typ `ch5ex*_online_query.m`, welche einen zufälligen Parameter erzeugt, sowohl die Truth- als auch die RB-Lösung bestimmt und anschließend tatsächlichen Fehler und a posteriori-Fehlerschätzer ausgibt.

### Beispiel 1

Um die Abbildungen und Tabellen nachzustellen, reicht es aus, das Skript `ch5ex1_plot_everything.m` auszuführen. 
Dieses lädt die bereits vorberechneten Daten und erzeugt die Abbildungen.

Will man die Daten neu berechnen, dann sollte dafür die Datei `ch5ex1_precompute_example.m` verwendet werden.
Diese führt zunächst `ch5ex1_run.m` aus, was die Offline-Phasen der SCM und RBM ausführt.
Anschließend wird `ch5ex1_run_analysis.m` ausgeführt, was weitere Berechnungen durchführt die zur Erzeugung der Abbildungen notwendig sind.
Die letzte verbleibende Datei `ch5ex1_problem_data.m` enthält die Modelldaten des Beispiels.

### Beispiel 2

Auch hier können mittels `ch5ex2_plot_everything.m` die Abbildung und die Tabelle nachgestellt werden.

Um die Berechnungen komplett neu auszuführen, sollte `ch5ex2_precompute_example.m` verwendet werden, welche `ch5ex2_run.m` ausführt und die generierten Daten speichert.
Erneut finden sich in `ch5ex2_problem_data.m` die Modelldaten des Beispiels.
