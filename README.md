# chrysper

Generalized Linear Model (family Negative Binomial)

Regression: 
`geneKO read counts ~ treatment*time + control variables + offset + ε`

A treatment that kills the tumor gradually is better for the patient than a treatment that kills the tumor too fast.

In this line, we want to detect with which gene knockouts there is a strong interaction between a (drug) treatment and its duration -- and if so, check the changes in cell death rate (or positive selection) along time (steady or sudden changes), whether there is monotonicity or not, etc.


# Cat
Tenim unes dades de CRISPR/Cas9 screening, que es resumeixen en uns "counts" dels quals l'abundància es correspon amb quant "dolent" or "beneficiós" és perdre un determinat gen per les cèl.lules (tumorals) d'un cultiu. És a dir, si inactivar un determinat gen a una d'aquestes cèl.lules és dolent per ella i la seva descendència (perquè per exemple pot ser un gen important per a la seva supervivència), doncs aquesta cèl.lula morirà o tindrà poca descendència, i llavors als resultats de l'anàlisi veurem pocs "counts" procedents d'aquesta cèl.lula (perquè n'hi haurà poques descendents al cultiu) comparat amb els counts de la resta de cèl.lules que tenen altres gens inactivats.

A més, tenim dades de cultius els quals s'han tractat amb un fàrmac antitumoral, i d'altres sense tractament (control).

Finalment, tenim mostres de cada cultiu a temps successius, així podem veure com les proporcions de counts van canviant al llarg del temps (potser la perdua d'un gen triga en ser dolenta, aleshores només quan han passat força dies podem detectar un descens en la quantitat de counts que representen les cèl.lules en les quals s'ha inactivat aquest gen).

Les variables de la NB regression poden ser categòriques, com la presència o no d'un fàrmac al cultiu, el tipus cel.lular si n'hi hagués més d'u... i time és el dia de mostreig (dia1, dia5, dia10...). Els asteriscs indiquen que volem testar les interaccions entre les variables.
Actualment aquesta regressió segueix un "Generalized Linear Model", de tipus "Negative Binomial" (com Poisson, perquè son counts, però amb més variància), que en R es pot indicar amb la funció 'glm.nb'. Crec que també s'estava provant amb lmer, que permeteix introduir mixed effects.

Així, la cosa seria que amb aquesta anàlisi es podria veure si, per exemple, hi ha una interacció entre un tractament antitumoral i la durada del tractament. Si és així, veuriem que la disminució de cèl.lules tumorals quan hi ha presència del tractament (menys counts totals, comparat amb el control) depèn de la durada del mateix. També podriem veure com de ràpid actua i altres aspectes.
