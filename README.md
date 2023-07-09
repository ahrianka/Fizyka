**OPIS REPOZYTORIUM**

W danym repozytorium jest umieszczony plik z kodem w języku Python dla uruchomienia symulacji problemów dwóch oraz trzech ciał.

Najpierw definiujemy stałe oraz wektory: stała grawitacyjna, masy ciał(gwiazd), położenia początkowe, prędkości początkowe i td. Potem z pomocą funkcji _**TwoBodyEquations**_ oraz _**ThreeBodyEquations**_ opisujemy równania ruchu systemu dwóch i trzech ciał odpowiednio. Zakładamy, że na ciała działają wyłącznie siły grawitacji innych ciał.

Funkcja _**update_animation**_ na podstawie równań ruchu z dwóch poprzednich funkcji rysuje drogę wszystkich ciał, którą pokonali od samego początku do pewnego momentu w czasie.

Potem z pomocą komendy _**plt.show()**_ jest wywoływana symulacja.
