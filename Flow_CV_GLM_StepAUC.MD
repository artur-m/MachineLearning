### Training Flow with Cross-Validation using GLM with categorization but without regularization 

I try to describe my feeling how training flow should looks like.

First In Polish.... ( send message if you are interested in English version )

## Założenia i Cel

Celem jest zbudowanie modelu scoringowego opartego na GLM. 
Dane wejściowe stanowi ramka danych (dataset):
  - w której jedna kolumna oznaczona jest jako target i przyjmuje wartości ze zbioru $$ \{0,1\}$$.
  - pozostałe są cechami kategorycznymi lub liczbowymi, które można będzie wykorzystać w produkcyjnym dzianianiu modelu z możliwością wystąpienia braków danych.
  
Zakładamy, że brak danych jest taką samą informacją jak inne podane wartości.
Wykonana zostanie kategoryzacja z transformacją. Przez tą operację mam na myśli ewentualne złączenie kategorii/ustalenie przedziałów na cechach ciągłych a następnie zakodowanie ramki danych do postaci macierzy liczbowej. Innymi słowy 2 ważne operacje typowa kategoryzacja + zakodowanie. Kodowaniem może by przejście na dummy variables lub WoE lub współczynniki w jednowymiarowym modelu lub jeszcze inne metody. Co warto odnotować zakładamy tutaj możliwość skorzystania z targetu przy poszukiwaniu optymalnej kategoryzacji. Założmy w tym miejscu, że procesem kategoryzacji możemy sterować poprzez wartości z przestrzeni $C$.

Ponieważ chcemy zastosować model liniowy ( bez żadnej regularyzacji ) więc przestrzenią parametrów definicujacą możliwe do zbudowania modele jest zbiór wszystkich podzbiorów cech z jakich skorzystamy. 
Stajemy więc przed optymalizacją w postaci optymalnej kategoryzacji + feature selection.

W tym miejscu chciałbym przedstawić propozycję FLOW poszukiwania i budowy modelu, który będzie miał najlepsze parmatry wg wybranej statystyki (AUC/lift/logloss) a jednocześnie zachowa zdolności predykcyjne. Operać się będzie o Cross-Validation + Test Set.

1. Podzielić dataset na dwie części: treningową(80%) i testową(20%)
2. Wybrać k ( propozycja 5 lub 10) i losowo podzielić zbiór treningowy na k części
3. Określić wstępny zbiór parametrów $CC \in C$, określających kategoryzację.
4. Dla każdego $c \in CC$
    - Wykonać niezależne kategoryzację wszystkich części wyszczególnionych w kroku 2. To jest przyjmując do nauki zbiór składający się z k-1 częsci, a ostatnią tylko do testowania.
    - Przechodzimy do poszukiwania cech. Przykładem może być tutaj Forward Step wg wybranej statystyki, jednak nic nie stoii na przeszkodzie by zastosować bardziej złożony algorytm przeszukujący szerszą przestrzeń podzbiorów cech. 
    - W każdym kroku kiedy testujemy dany układ cech budując k modeli i  wyznaczając rozkład statystyki na podstawie k testów na fragmentach pozostawionych do testów.
    - Ważnym krokiem jest sposób w jaki podejmiemy decyzję, który model jest lepszy na podstawie dwóch rozkładów statystyki. Tym co chcielibyśmy uzyskać to jaknajkorzystniejsza wartość średnia i jednocześnie jak najmnniejszy rozrzut wartości. Duży rozrzut świadczy o przeuczaniu się modelu.
5. Po zakończeniu tych zautomatyzowanych poszukiwań powinniśmy dostać zbiór cech i kategoryzację na nich, która prowadzi do stabilnego modelu. Możemy wtedy na całym zbiorze treningowym zbudować model z wyznaczonymi parametrami i przetestować na zbiorze testowym. Oczekujemy zobaczym zbliżoną wartość statystyki jak uśrtednione podczas CV. Jeśli tak się nie stanie należy zdecydować czy może należy rozpocząć budowę od początku z wykorzystaniem innych parametrów/warunków stopu/oceny który model lepszy czy może zaufać wynikom CV i liczyć, że złe działanie na zbiorze testowym to efekt małej próby.
