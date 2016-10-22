


#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
#---------------------------------------------------------------------- Historia zmian -------------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------#


## Nr wersji:          04
## Zakres zmian:
# Umo?liwienie obs?ugi targetu z przedzia?u [0,1] (uog?lnienie)
#############################################################################################################################################################################


## Nr wersji:          05
## Zakres zmian:
# Dodanie nowej miary jako?ci podzia?u zmiennej - zmienno?? wariancji ("deviance")
#############################################################################################################################################################################


## Nr wersji:          06
## Zakres zmian:
# Dodanie nowej miary jako?ci podzia?u zmiennej - zmienno?? parametru phi w regresji beta ("phi_vol")
# Korekta miary "deviance": globalne odchylenie sum(c_c*c_sd)/sum(c_c) zast?pujemy przez sd_total = sd(target). 
#############################################################################################################################################################################




#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
#---------------------------------------------------------- SPECYFIKACJA FUNKCJI 'var_partition_fitting_measure' -----------------------------------------------------------#
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------#


# ----- CHARAKTERYSTYKA FUNKCJI --------------------------------------------------------------------------------------------------------------------------------------------#
# Funkcja oblicza miar? dyskryminacji dla zmiennej wejsciowej. Obecnie dopuszczalne s? dwie miary: IV lub p-value testu chi2. Funkcja umozliwia obs?ug? wyj?tkowych przypadk?w,
# kt?re komplikuj? obliczenie IV poprzez zastosowanie tzw sztucznej licznosci dla przedzia?u, kt?rego rzeczywista licznos? jest r?wna 0.
# FUNKCJA JEST UOG?LNIENIEM var_partition_fitting_measure (target z przedzia?u [0,1] a wi?c mozliwe s? u?amkowe liczno?ci; st?d FC -fractional counts)


# ----- ARGUMENTY WEJSCIOWE ------------------------------------------------------------------------------------------------------------------------------------------------#
# variable         -> zmienna dla kt?rej wyliczana jest miara
# target           -> target - rzeczywiste wykonanie - ZMIENNA PRZYJMUJ?CA WARTO?CI Z PRZEDZIA?U [0,1]. Interpretacja: 0 -->  w 100% "good", w 0% "bad";
#                     0.4 --> w 40% "good", w 60% '"bad"; 1 --> w 0% "good", w 100% '"bad"
# flag_artif_count -> warto? logiczna (TRUE/FALSE) informuj?ca czy w przypadku problem?w z obliczeniem IV stosowa? sztuczn? licznos?.
# artificial_count -> wartos? dla sztucznej licznosci



# ----- OBIEKTY ZWRACANE ---------------------------------------------------------------------------------------------------------------------------------------------------#
# var_partition_fitting_measure -> obliczona miara




#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------#



var_partition_fitting_measure  = function( variable, target, measure = "iv", flag_artif_count = FALSE, artificial_count = 0.5)
{
  
  ####tab = table(variable, target)
  tab = as.table( cbind(tapply(1-target, variable, FUN = "sum"), tapply(target, variable, FUN = "sum")) )
  colnames(tab) = c(0,1)
  tab[which(is.na(tab))] = 0 
  if( flag_artif_count == TRUE )
  { tab[which( tab == 0 )] = artificial_count }
  if( measure == "iv" )
  {
    var_tab_prob        = prop.table(tab, margin=2)
    woe                 = as.numeric(log( var_tab_prob[,1] / var_tab_prob[,2] ))
    
    ind_error           = which(is.element(woe, c(Inf,-Inf, NA, NaN)))
    if( length(ind_error) > 0 )
    { stop("Error during IV measure calculation!") }
    var_partition_fitting_measure = sum( (var_tab_prob[,1] - var_tab_prob[,2])*woe )
  }
  else if( measure == "chisq" )
  { var_partition_fitting_measure = chisq.test(tab)$p.value }
  else if( measure == "deviance" )  # zmienno?? wariancji - istotne dla modelowania zmienno?ci w regresji BETA
  {
    c_sd                          = tapply(target,variable,sd)
    c_sd[which(is.na(c_sd) | c_sd == Inf )] = sd(target)
    sd_total                      = sd(target)
    c_c                           = tapply(target,variable,length)
    var_partition_fitting_measure = sqrt(sum(c_c*(c_sd-sd_total)^2)/sum(c_c))
  }
  else if( measure == "phi_vol")  # volatility of phi
  {
    phi_total                     = mean(target)*(1-mean(target))/var(target) - 1
    phi                           = tapply(target,variable,mean)*( 1 - tapply(target,variable,mean) )/tapply(target,variable,var) - 1
    phi[which(is.na(phi) | phi == Inf )] = phi_total   
    p                             = prop.table(table(variable))
    var_partition_fitting_measure = sum(p*(phi-phi_total)^2)
  }
  else
  { stop("Undefined measure type!") }
  return( var_partition_fitting_measure )
}





#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
#---------------------------------------------------------------------- Historia zmian -------------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------#


## Nr wersji:          11
## Zakres zmian:
# 1. zmiana nazwy zapisywanego pliku (skr?cenie)
#     OBIEKT          STARA NAZWA              NOWA NAZWA
#     plik csv        __partition_process      __pproc
#############################################################################################################################################################################


## Nr wersji:          12
## Zakres zmian:
# Umo?liwienie obs?ugi targetu z przedzia?u [0,1] (uog?lnienie)
#############################################################################################################################################################################


## Nr wersji:          13
## Zakres zmian:
# Korekta b??du w obs?udze parametru 'constrains' dla zmiennych dyskretnych, przy podaniu kilku wektor?w dochodzi?o do pomieszania etykiet
#############################################################################################################################################################################


## Nr wersji:          14
## Zakres zmian:
# Korekta b??du w obs?udze parametru 'constrains' dla zmiennych dyskretnych, przy podaniu klas jedno-elementowych p?tla ??cz?ca etykiety generowa?a b??d
#############################################################################################################################################################################



#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
#---------------------------------------------------------- SPECYFIKACJA FUNKCJI 'var_partition' ---------------------------------------------------------------------------#
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------#


# ----- CHARAKTERYSTYKA FUNKCJI --------------------------------------------------------------------------------------------------------------------------------------------#
# Funkcja dokonuje partycji zmiennej wykorzystuj?c prosty algorytm optymalizacji zadanej miary mocy dyskryminacyjnej. Algorytm zaczyna dzia?anie od okreslonej przez uzytkownika
# liczby podzia??w lub wprost zadanych podzia??w. W przypadku zmiennej dyskretnej na wejsciu algorytmu mamy wszystkie wyst?puj?ce kategorie zmiennej, chyba ze wykorzystano opcj?
# w postaci argumentu 'constrains', kt?ry umozliwia scalenie okreslonych kategorii juz na pocz?tku dzia?ania algorytmu.
# W przypadku zmiennej ci?g?ej gdy uzytkownik wprowadza liczb? podzia??w (nie zas konkretny podzia?) funkcja nie gwarantuje unikni?cia problem?w z zerow? licznosci? uzyskanego
# podzia?u w kt?rejs z podpopulacji wyznaczonej przez argument 'target'. Gdy podzia? jest uzytkownika problem ten jest rozwiazywany przez "sztuczn? licznos? (argument
# 'artificial_count'). W przypadku zmiennej dyskretnej problem ten jest rozwi?zany przez scalanie odpowiednich kategorii (patrz funkcja 'bands_merging').
# Spos?b dzia?ania algorytmu:
# 1. Algorytm dla zmiennej ci?g?ej polega na ??czeniu s?siaduj?cych przedzia??w (wyj?tkowo traktowany jest przypadek braku danych. Ta specyficzna grupa nie podlega ??czeniu z
#    zadnym przedzia?em). Dla zmiennej dyskretnej nast?puje przesortowanie polegaj?ce na ustawieniu kategorii wed?ug gradacji 'bad rates' (malej?co). Dla zmiennej dyskretnej
#    kategoria braku danych traktowana jest jako jedna z kategorii i nie podlega zadnym wyj?tkom.
# 2. Przy zadanym podziale algorytm sprawdza wszystkie pary s?siaduj?cych przedzia??w/kategorii (dla zmiennych dyskretnych 's?siadowanie' wyznacza porz?dek 'bads rate') i scala
#    w jeden przedzia?/kategori? t? par? kt?ra powoduje najmniejsz? strat? w stosowanej mierze dyskryminacji.
# 3. Algorytm zatrzymuje si? jesli wartos? miary zmieni si? o wi?cej niz zadana wartos? (argument 'eps') w stosunku do wyjsciowej wartosci miary lub gdy nie ma juz co scalac.
#
# Istnieje mozliwos? niestosowania optymalizacji (argument 'flag_optim'), w?wczas funkcja wykorzystywana jest do wprowadzenia nietypowych w?asnych podzia??w, na przyk?ad
# zawieraj?cych pojedynczy punkt, czy tez scalanie "izolowanych" przedzia??w (np. (0,1] i (3,7] - argument 'constrain' ale tylko dla zmiennej ci?g?ej!) lub tez unikniecia
# problem?w z funkcj? 'cut' gdy argument 'breaks' nie spe?nia za?oze? (cz?sty problem "breaks are not unique").




# ----- ARGUMENTY WEJSCIOWE ------------------------------------------------------------------------------------------------------------------------------------------------#
# main_path  -> sciezka do katalogu w kt?rym znajduj? si? podkatalogi zawieraj?ce poszczeg?lne zmienne do zakodowania. Podkatalogi to wynik funkcji '1D Analysis'
# data_set   -> ramka danych zawieraj?ca zmienne do zakodowania
# var_names  -> wektor z nazwami zmiennych, kt?re maja zostac zakodowane (opcja domyslna: kodowane s? wszystkie zmienne z 'data_set'

# variable         -> zmienna, kt?ra ma podlega? partycji
# flag_contin      -> warts? logiczna (TRUE/FALSE) informuj?ca czy zmienna jest ci?g?a
# var.name         -> nazwa zmiennej 'variable' (wykorzystywana przy zapisie wynik?w partycji)
# target,          -> target - rzeczywiste wykonanie - dziel?ce zmienn? 'variable' na dwie grupy (0-1)
# no_initial_part  -> liczba pocz?tkowych podzia??w
# user_breaks      -> wektor punkt?w zadaj?cych konkretny podzia? zmiennej 'variable'
# flag_user_breaks -> warts? logiczna (TRUE/FALSE) informuj?ca czy uzytkownik wprowadza w?asny podzia? (za pomoc? argumentu 'user_breaks'). Wprowadzenie tego argumentu wynika z
#                     koniecznosci rozr?znienia sytuacji, gdy pusty wektor 'user_breaks' oznacza, ze uzytkownik nie chce w?asnych podzia??w od tej gdy chce jeden przedzia? dla
#                     zmiennej (-Inf; Inf) (zwykle wiedz?c, ze przyb?dzie jeszcze jeden w postaci braku danych)
# eps              -> liczba okreslaj?ca jaki % miary dyskryminacji moze zosta? utracony na rzecz scalania przedzia??w/kategorii
# constrains       -> lista zawieraj?ca wektory przedzia??w/kategorii, kt?re maj? by? scalone. Kazdy wektor z listy stanowi ograniczenie wymuszaj?ce utworzenie z  kategorii
#                     zestawionych w wektorze jednej "duzej" kategorii kt?ra pojawi si? na wejsciu algorytmu optymalizacji. Dla zmiennej ci?g?ej 'constrains' rozumiany jest jako
#                     lista zawieraj?ca wektory przedzia??w, kt?re maj? by? scalone mimo iz nie s?siaduj? ze sob?. Taka opcja jest przewidziana tylko w sytuacji, gdy NIE
#                     zak?ada si? optymalizacji! (flag_optim = FALSE)
# threshold        -> pr?g licznosci ponizej kt?rego ma zosta? wykorzystana "sztuczna licznos?"
# measure_type     -> typ miary wykorzystany w algorytmie optymalizacji (obecnie IV lub p-value testu chi2)
# optimization     -> rodzaj optymalizacji (minimalizacja lub maksymalizacja) dla miary
# rounding         -> liczba miejsc po przecinku do kt?rej zaokr?glone zostan? punkty podzia?u. Opcja domyslna to brak zaokr?glenia
# flag_artif_count -> wartos? logiczna (TRUE/FALSE) informuj?ca czy ma by? wykorzystana "sztuczna licznos?"
# artificial_count -> wartos? ("sztuczna licznosc"), kt?r? ma by? zast?piona licznosc przedzia?u/kategorii kt?ra jest ponizej wartosci progowej ('threshold')
# flag_optim       -> wartos? logiczna (TRUE/FALSE) informuj?ca czy ma by? zastosowany algorytm optymalizacji
# write_flag       -> wartos? logiczna (TRUE/FALSE) informuj?ca czy wynik dzia?ania funkcji ma by? zapisany
# write_file_path  -> sciezka zapisu danych




# ----- OBIEKTY ZWRACANE ---------------------------------------------------------------------------------------------------------------------------------------------------#
#output - lista nast?puj?cych obiekt?w:
# variable       -> faktor, kt?ry zawiera uzyskany ko?cowy podzia? zmiennej variable_cut,
# bands          -> karta zawieraj?ca kod przewidziany do zakodowania zmiennej w spos?b identyczny kt?ry wygenerowa?a funkcja,
# bands_count    -> table kt?ry zawiera licznosci w poszczeg?lnych przedzia?ach/kategoriach
# measure_value  -> wartos? miary dyskryminacyjnej uzyskanej przy "optymalnym" podziale (jesli nie stosowano optymalizacji miara r?wniez jest zwracana)
# merging_raport -> raport podsumowuj?cy kolejne kroki w procesie optymalizacji



#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------#


var_partition = function(
  variable,
  flag_contin     = TRUE,
  var.name        = "variable",
  target,
  no_initial_part = 5,
  user_breaks     = c(),
  flag_user_breaks = FALSE,
  eps             = 0.1,
  constrains      = list(),
  threshold       = 0,
  measure_type    = "iv",
  optimization    = "max",
  rounding        = c(),
  flag_artif_count = FALSE,
  artificial_count = 0.5,
  flag_optim      = TRUE,
  write_flag      = TRUE,
  write_file_path = "d:\\documents and settings\\wosarn\\Moje dokumenty\\Workspace\\"
)
{
  # -- ustawienie flag
  flag_na       = FALSE
  flag_na_thres = FALSE
  flag_loop_end = FALSE
  
  #======================================================================================================================================#
  #============================= CONTINUOUS VARIABLE =======================================================================================#
  #======================================================================================================================================#
  if(flag_contin == TRUE)
  {
    ind_na      = which(is.na(variable))
    flag_na     = ifelse( length(ind_na) > 0, TRUE, FALSE)
    # -- sprawdzenie czy NA spe?nia warunek "threshold"
    if( flag_na == TRUE & flag_artif_count == FALSE )
    {
      variable_tmp          = vector( mode = "numeric", length = length(variable) )
      variable_tmp[-ind_na] = 1
      ####tab_tmp               = table(variable_tmp, target)
      tab_tmp = as.table( cbind(tapply(1-target, variable_tmp, FUN = "sum"), tapply(target, variable_tmp, FUN = "sum")) ) #uog?lnienie
      colnames(tab_tmp) = c(0,1)
      tab_tmp[which(is.na(tab_tmp))] = 0  
      if( length( which(tab_tmp[1,] <= threshold) ) > 0 )
      {
        flag_na_thres = TRUE
        warning("Missing values are not taken into account during optimization process" )
      }
      rm(variable_tmp)
      rm(tab_tmp)
    }
    loop_ind = no_initial_part+1
    # -- ustalanie liczby podzia??w, tak aby spe?niony by? warunek "threshold"
    while( flag_loop_end == FALSE )
    {
      if( flag_user_breaks == FALSE )
      {
        q                = as.numeric( quantile( variable, seq(0, 1, length=loop_ind), na.rm = TRUE ) )
        if( length(rounding) == 1 & length(q) > 2)
        { q[2:(length(q)-1)] = round(q[2:(length(q)-1)], digits = rounding ) }
      }
      else
      {
        if( length(user_breaks) > 0 )
        {
          if( length(rounding) == 1)
          {
            ind_not_extr = which( !is.element(user_breaks, c(min(variable, na.rm = TRUE),max(variable, na.rm = TRUE))) )
            if( length(ind_not_extr)>0 )
            { user_breaks[ind_not_extr] = round(user_breaks[ind_not_extr], digits = rounding ) }
          }
          ind_atoms = which(table(user_breaks) > 1)
          atoms     = as.numeric( names(table(user_breaks))[ind_atoms] )
        }
        else
        {
          ind_atoms = integer(0)
          atoms     = numeric(0)
        }
        q_aux = unique( c(min(variable, na.rm = TRUE), user_breaks, max(variable, na.rm = TRUE)) )
        q     = sort(c(q_aux, atoms))
        #        q                = user_breaks
        flag_artif_count = TRUE
      }
      q_tmp                       = sort(unique(q))
      tq                          = table(q)
      ind_multip                  = which( as.numeric(tq)>1 )
      var_cut_tmp                 = cut(variable, breaks = q_tmp, right = TRUE, include.lowest = TRUE, dig.lab = 6)
      bands_number_aux            = length(q_tmp)-1
      bands                       = as.data.frame( matrix( 0, nrow = bands_number_aux, ncol = 4) )
      names(bands)                = c("low_cl", "low_val", "upp_val", "upp_cl") #, "label")
      #     bands$labels                = ""
      bands$labels                = NA
      bands[1:bands_number_aux,2] = q_tmp[1:(length(q_tmp)-1)]
      bands[1:bands_number_aux,3] = q_tmp[2:length(q_tmp)]
      bands[1,c(1,4)]             = 1
      if( bands_number_aux > 1 )
      {
        bands[2:bands_number_aux,4] = 1
        bands[2:bands_number_aux,1] = 0
      }
      if( length( ind_multip ) == 0 )
      { variable_cut = factor( var_cut_tmp, levels = levels(var_cut_tmp), exclude = NULL) }
      # { variable_cut = factor( var_cut_tmp, exclude = NULL) }
      else        # szczeg?lne traktowanie atom?w, skupiaj?cych zbyt duz? mas? rozk?adu (takie atomy utworz? osobne kategorie)
      {
        levels_aux = levels(var_cut_tmp)
        for( i in 1 : length(ind_multip) )
        {
          correction = ind_multip[i] + i-1
          if( ind_multip[i] == 1 )
          {
            levels_aux                 = c( q_tmp[ind_multip[i]],levels_aux)
            bands                      = as.data.frame(rbind(c(1,q_tmp[ind_multip[i]],q_tmp[ind_multip[i]],1, NA),bands))
            bands[2,1]                 = 0
          }
          #  else if( ind_multip[i] == length(ind_multip) )
          else if( ind_multip[i] == length(q_tmp) )
          {
            levels_aux                 = c( levels_aux, q_tmp[ind_multip[i]])
            bands                      = as.data.frame(rbind(bands,c(1,q_tmp[ind_multip[i]],q_tmp[ind_multip[i]],1, NA)))
            bands[(dim(bands)[1]-1),4] = 0
          }
          else
          {
            levels_aux                 = c( levels_aux[1:(correction-1)], q_tmp[ind_multip[i]], levels_aux[correction:length(levels_aux)] )
            bands                      = as.data.frame(rbind(bands[1:(correction-1),], c(1,q_tmp[ind_multip[i]],q_tmp[ind_multip[i]],1, NA),bands[correction:dim(bands)[1],]))
            bands[(correction-1),4]    = 0
          }
        }
        variable_cut = factor( var_cut_tmp, levels = levels_aux, exclude = NULL)
        for( j in 1 : length(ind_multip) )
        {
          ind               = which( variable == q_tmp[ind_multip[j]] )
          variable_cut[ind] = q_tmp[ind_multip[j]]
        }
      }
      if( flag_na == TRUE )
      {
        variable_cut                     = factor( variable_cut, levels = c(levels(variable_cut),NA), exclude = NULL)
        ind_lev_na                       = which(is.na(levels(variable_cut)))
        levels(variable_cut)[ind_lev_na] = "blank"
        bands                            = as.data.frame(rbind(c(1, NA, NA, 1, NA),bands))
        variable_cut                     = relevel(variable_cut,"blank")
      }
      ####tab  = table(variable_cut, target)
      tab = as.table( cbind(tapply(1-target, variable_cut, FUN = "sum"), tapply(target, variable_cut, FUN = "sum")) ) #uog?lnienie
      colnames(tab) = c(0,1)
      tab[which(is.na(tab))] = 0      
      if( flag_na_thres == TRUE )
      { tab = tab[-1,] }
      if( length( which(tab <= threshold) ) > 0 & flag_artif_count == FALSE )
      { loop_ind = loop_ind - 1 }
      else
      { flag_loop_end = TRUE  }
      if( loop_ind == 1 )      ## czy nie loop_ind == 1 lub nawet loop_ind == 2 ???
      { stop("Error during merging process") }
      labels_merging = levels(variable_cut)
    }
    if( ( loop_ind < (no_initial_part+1) | length( ind_multip ) > 0 ) & flag_user_breaks == FALSE)
    { warning("Initial partition has been changed in authomatic optimization process") }
    rm(loop_ind)
    
    # -- szukanie optymalnego podzialu
    ###  initial_bands_number = length( labels_merging )
    bands_number = length(labels_merging)
    initial_ind  = 1
    if( flag_na == TRUE )
    { initial_ind = 2 }
    
    loop_ind                    = 1
    merging_summary             = as.data.frame( rbind(c(1:4,1:bands_number+4)) )
    names(merging_summary)      = c("bands_number", "measure_type", "measure_value", "epsilon", paste("band", 1:bands_number, sep = "_") )
    merging_summary[loop_ind,1] = bands_number
    merging_summary[loop_ind,2] = measure_type
    if( flag_na_thres == TRUE & flag_artif_count == FALSE)
    {
      variable_cut_tmp            = factor( variable_cut[-which(variable_cut == "blank")] )
      target_tmp                  = target[-which(variable_cut == "blank")]
      merging_summary[loop_ind,3] = var_partition_fitting_measure( variable = variable_cut_tmp, target = target_tmp, measure = measure_type ) #uog?lnienie
    }
    else
    { merging_summary[loop_ind,3] = var_partition_fitting_measure( variable = variable_cut, target = target, measure = measure_type, flag_artif_count = flag_artif_count, artificial_count = artificial_count)  } #uog?lnienie
    merging_summary[loop_ind,4]                = eps
    merging_summary[loop_ind,1:bands_number+4] = labels_merging
    
    measure_comparision_value = merging_summary[loop_ind,3]
    loop_ind                  = loop_ind + 1
    #    flag_loop_end             = ifelse( initial_ind == bands_number, TRUE, FALSE)
    flag_loop_end             = ifelse( flag_optim == FALSE, TRUE, ifelse( initial_ind == bands_number, TRUE, FALSE) )
    
    while( flag_loop_end == FALSE  )
    {
      measure_tmp  = vector("numeric", length = bands_number-1) * NA
      for( j in initial_ind : (bands_number-1) )
      {
        variable_cut_tmp                   = variable_cut
        label_new                          = paste( labels_merging[j], labels_merging[j+1], sep = "_")
        ind_lev1                           = which( levels(variable_cut_tmp) == labels_merging[j] )
        ind_lev2                           = which( levels(variable_cut_tmp) == labels_merging[j+1] )
        levels(variable_cut_tmp)[ind_lev1] = label_new
        levels(variable_cut_tmp)[ind_lev2] = label_new
        if( flag_na_thres == TRUE & flag_artif_count == FALSE)
        {
          variable_cut_tmp2 = factor( variable_cut_tmp[-which(variable_cut_tmp == "blank")] )
          target_tmp        = target[-which(variable_cut_tmp == "blank")]
          measure_tmp[j]    = var_partition_fitting_measure( variable = variable_cut_tmp2, target = target_tmp, measure_type) #uog?lnienie
        }
        else
        { measure_tmp[j]    = var_partition_fitting_measure( variable_cut_tmp, target, measure_type, flag_artif_count = flag_artif_count, artificial_count = artificial_count) } #uog?lnienie
      }
      if( optimization == "max" )
      {
        ind_opt           = which.max( measure_tmp )
        measure_condition = (measure_comparision_value - measure_tmp[ind_opt]) / measure_comparision_value
      }
      else
      {
        ind_opt           = which.min( measure_tmp )
        measure_condition = (measure_tmp[ind_opt] - measure_comparision_value) / measure_comparision_value
      }
      if( measure_condition <= eps )
      {
        ind_lev1_opt                                         = which( levels(variable_cut) == labels_merging[ind_opt] )
        ind_lev2_opt                                         = which( levels(variable_cut) == labels_merging[ind_opt+1] )
        levels(variable_cut)[ind_lev1_opt]                   = paste( labels_merging[ind_opt], labels_merging[ind_opt+1], sep = "_")
        levels(variable_cut)[ind_lev2_opt]                   = paste( labels_merging[ind_opt], labels_merging[ind_opt+1], sep = "_")
        labels_merging                                       = levels(variable_cut)
        merging_summary[loop_ind,1]                          = length(labels_merging)
        merging_summary[loop_ind,2]                          = measure_type
        merging_summary[loop_ind,3]                          = measure_tmp[ind_opt]
        merging_summary[loop_ind,4]                          = eps
        merging_summary[loop_ind,1:length(labels_merging)+4] = labels_merging
        bands[ind_lev1_opt, 3:4]                             = bands[ind_lev2_opt, 3:4]
        bands                                                = bands[-ind_lev2_opt,]
      }
      loop_ind     = loop_ind + 1
      bands_number = length(labels_merging)
      if( initial_ind == bands_number | measure_condition > eps )
      { flag_loop_end = TRUE }
      rm(measure_tmp)
    }
    
    bands[initial_ind,1:2] = c(0,NA)
    bands[length(labels_merging),3:4] = c(NA,0)
    brackets = c( "(", "[", ")", "]" )
    for( i in 1:dim(bands)[1] )
    {
      if( is.na(bands[i,2]) & is.na(bands[i,3]) & sum(bands[i,c(1,4)]) == 2 )
      { bands$labels[i] = "blank" }
      else if(  is.na(bands[i,2]) & is.na(bands[i,3]) & sum(bands[i,c(1,4)]) == 0 )
      { bands$labels[i] = "(-Inf; Inf)"}
      else if( is.na(bands[i,2]) )
      { bands$labels[i] = paste( paste( paste( paste( brackets[bands[i,1]+1], -Inf, sep="" ),  ";", sep="" ), bands[i,3], sep=" " ), brackets[bands[i,4]+3], sep="" ) }
      else if( is.na(bands[i,3]) )
      { bands$labels[i] = paste( paste( paste( paste( brackets[bands[i,1]+1], bands[i,2], sep="" ),  ";", sep="" ), Inf, sep=" " ), brackets[bands[i,4]+3], sep="" ) }
      else if( bands[i,2] == bands[i,3] )
      { bands$labels[i] = paste(bands[i,2]) }
      else
      { bands$labels[i] = paste( paste( paste( paste( brackets[bands[i,1]+1], bands[i,2], sep="" ),  ";", sep="" ), bands[i,3], sep=" " ), brackets[bands[i,4]+3], sep="" )  }
    }
    levels(variable_cut) = bands$labels
    
    if(length(constrains)>0 & flag_optim == FALSE)
    {
      for( i in 1: length(constrains) )
      {
        #  ind_constr = which( is.element(part_card[,5],constrains[[i]]) )
        ind_constr1 = which( is.element(bands$labels,constrains[[i]]) )
        ind_constr2 = which( is.element(levels(variable_cut),constrains[[i]]) )
        if( length(ind_constr1) == length(constrains[[i]]) & length(constrains[[i]]) >1 & length(ind_constr2) == length(constrains[[i]]) )
        {
          new_label                         = paste(bands$labels[ind_constr1], collapse="_")
          bands$labels[ind_constr1]         = new_label
          levels(variable_cut)[ind_constr2] = new_label
        }else
        {warning(paste("problems with constrains:", paste(constrains[[i]],collapse=",")))}
      }
    }
    
    
    
    
  }
  #======================================================================================================================================#
  #============================= DISCRET VARIABLE =======================================================================================#
  #======================================================================================================================================#
  else
  {
    variable_cut                 = factor(variable, exclude = NULL)
    ind_na                       = which(is.na(levels(variable_cut)))
    levels(variable_cut)[ind_na] = "blank"
    tapp                         = tapply(target, variable_cut, FUN = "mean")
    ind_order                    = order( as.numeric(tapp) )
    labels_merging               = levels(variable_cut)[ind_order]
    
    bands                        = as.data.frame( matrix( NA, nrow = length(labels_merging), ncol = length(labels_merging)+1) )
    bands[,1]                    = labels_merging
    names(bands)                 = c(paste("atr",1:length(labels_merging),sep=""),"label") #labels_merging
    bands$label                  = labels_merging
    
    if( length(constrains) > 0 )
    {
      rows_to_remove=integer(0)
      for( i in 1 : length(constrains) )
      {
        rep1                              = matrix(rep( labels_merging        , times = length(constrains[[i]])), nrow = length(labels_merging), ncol = length(constrains[[i]]))
        rep3                              = matrix(rep( paste(constrains[[i]]), each  = length(labels_merging)) , nrow = length(labels_merging), ncol = length(constrains[[i]]))
        ind_constr1                       = which( rep1 == rep3, arr.ind = TRUE)[,1]
        rep2                              = matrix(rep( levels(variable_cut)  , times = length(constrains[[i]])), nrow = length(levels(variable_cut)), ncol = length(constrains[[i]]))
        rep3                              = matrix(rep( paste(constrains[[i]]), each  = length(levels(variable_cut))) , nrow = length(levels(variable_cut)), ncol = length(constrains[[i]]))
        ind_constr2                       = which( rep2 == rep3, arr.ind = TRUE)[,1]
        ##          bands[,ind_constr1[1]]            = constrains[[i]]
        bands[ind_constr1[1],1:(length(ind_constr1))] = labels_merging[ind_constr1]
        label_new                         = paste(labels_merging[ind_constr1],collapse="_")
        bands$label[ind_constr1[1]]       = label_new
        ###bands                             = bands[-ind_constr1[2:length(ind_constr1)],] ###usuni?cie zbednych wierszy wykonujemy po zako?czeniu p?tli
        if(length(ind_constr1) > 1) { #tylko w takim wypadku s? wiersze do usuni?cia
          rows_to_remove=c(rows_to_remove,ind_constr1[2:length(ind_constr1)])###zbieramy tylko indeksy wierszy do usuni?cia
        }
        #          names(bands)[ind_constr1[1]]      = label_new
        levels(variable_cut)[ind_constr2] = label_new
      }
      
      if(length(rows_to_remove) > 0) {
        bands <- bands[-rows_to_remove, ] ##usuni?cie zb?dnych wierszy (zosta?y po??czone z innymi)
      }
      
      cols_to_remove <- which(apply(bands, 2, function(x){sum(!is.na(x))}) == 0)
      if(length(cols_to_remove) >0) {
        bands <- bands[, ] ## usuni?cie zb?dnych kolumn (zawieraj? wy??cznie NA)
      }
      
      tapp           = tapply(target, variable_cut, FUN = "mean")
      ind_order      = order( as.numeric(tapp) )
      labels_merging = levels(variable_cut)[ind_order]
    }
    length_tmp                   = length( unique(variable_cut) )
    if( flag_artif_count == FALSE )
    {
      list_tmp       = bands_merging( variable = variable_cut, target = target, threshold = 0, labels_merged = bands)  #uog?lnienie
      variable_cut   = list_tmp$variable
      bands          = list_tmp$labels_merged
      tapp           = tapply(target, variable_cut, FUN = "mean")
      ind_order      = order( as.numeric(tapp) )
      labels_merging = levels(variable_cut)[ind_order]
      rm(list_tmp)
    }
    if( length_tmp > length(labels_merging) )
    { warning( "Initial partition in authomatic optimization process has been changed" ) }
    rm( length_tmp )
    
    ##-- szukanie optymalnego podzialu
    bands_number                = length(labels_merging)
    initial_ind                 = 1
    loop_ind                    = 1
    merging_summary             = as.data.frame( rbind(c(1:4,1:bands_number+4)) )
    names(merging_summary)      = c("bands_number", "measure_type", "measure_value", "epsilon", paste("band", 1:bands_number, sep = "_") )
    merging_summary[loop_ind,1] = bands_number
    merging_summary[loop_ind,2] = measure_type
    #merging_summary[loop_ind,3] = var_partition_fitting_measure( variable = variable_cut, target = target, measure = measure_type, flag_artif_count = flag_artif_count, artificial_count = artificial_count)
    merging_summary[loop_ind,3] = var_partition_fitting_measure( variable = variable_cut, target = target, measure = measure_type, flag_artif_count = flag_artif_count, artificial_count = artificial_count) #uog?lnienie
    merging_summary[loop_ind,4] = eps
    merging_summary[loop_ind,1:bands_number+4] = labels_merging
    
    measure_comparision_value = merging_summary[loop_ind,3]
    loop_ind                  = loop_ind + 1
    #    flag_loop_end             = ifelse( initial_ind == bands_number, TRUE, FALSE)
    flag_loop_end             = ifelse( flag_optim == FALSE, TRUE, ifelse( initial_ind == bands_number, TRUE, FALSE) )
    
    while( flag_loop_end == FALSE  )
    {
      measure_tmp  = vector("numeric", length = bands_number-1) * NA
      for( j in initial_ind : (bands_number-1) )
      {
        variable_cut_tmp                   = variable_cut
        label_new                          = paste( labels_merging[j], labels_merging[j+1], sep = "_")
        ind_lev1                           = which( levels(variable_cut_tmp) == labels_merging[j] )
        ind_lev2                           = which( levels(variable_cut_tmp) == labels_merging[j+1] )
        levels(variable_cut_tmp)[ind_lev1] = label_new
        levels(variable_cut_tmp)[ind_lev2] = label_new
        measure_tmp[j]    = var_partition_fitting_measure( variable_cut_tmp, target, measure_type, flag_artif_count = flag_artif_count, artificial_count = artificial_count) #uog?lnienie
      }
      if( optimization == "max" )
      { ind_opt = which.max( measure_tmp ) }
      else
      { ind_opt = which.min( measure_tmp ) }
      
      if( optimization == "max" )
      { measure_condition = (measure_comparision_value - measure_tmp[ind_opt]) / measure_comparision_value }
      else
      { measure_condition = (measure_tmp[ind_opt] - measure_comparision_value) / measure_comparision_value  }
      
      if( measure_condition <= eps )
      {
        ind_lev1_opt                                         = which( levels(variable_cut) == labels_merging[ind_opt] )
        ind_lev2_opt                                         = which( levels(variable_cut) == labels_merging[ind_opt+1] )
        levels(variable_cut)[ind_lev1_opt]                   = paste( labels_merging[ind_opt], labels_merging[ind_opt+1], sep = "_")
        levels(variable_cut)[ind_lev2_opt]                   = paste( labels_merging[ind_opt], labels_merging[ind_opt+1], sep = "_")
        ind_order                                            = order( as.numeric(tapply(target, variable_cut, FUN = "mean")) )
        labels_merging                                       = levels(variable_cut)[ind_order]
        merging_summary[loop_ind,1]                          = length(labels_merging)
        merging_summary[loop_ind,2]                          = measure_type
        merging_summary[loop_ind,3]                          = measure_tmp[ind_opt]
        merging_summary[loop_ind,4]                          = eps
        merging_summary[loop_ind,1:length(labels_merging)+4] = labels_merging
        ind1                                                 = which( is.na(bands[ind_opt,]) )[1]
        ind2                                                 = which( is.na(bands[ind_opt+1,]) )[1]
        bands[ind_opt, ind1:(ind1+ind2-2)]                   = bands[ind_opt+1, 1:(ind2-1)]
        bands$label[ind_opt]                                 = paste( bands$label[ind_opt], bands$label[ind_opt+1], sep="_")
        #        names(bands)[ind_opt]                                = paste( names(bands)[ind_opt], names(bands)[ind_opt+1], sep="_")
        bands                                                = bands[-(ind_opt+1),]
      }
      loop_ind     = loop_ind + 1
      bands_number = length(labels_merging)
      if( initial_ind == bands_number | measure_condition > eps )
      { flag_loop_end = TRUE }
      rm(measure_tmp)
    }
    ind_first_not_na = max(which( !is.na(bands[,-dim(bands)[2]]), arr.ind = TRUE)[,2])
    bands            = as.data.frame( bands[,c(1:ind_first_not_na, which(names(bands) == "label"))] )
    names(bands)     = c(paste("atr",1:(dim(bands)[2]-1),sep=""),"label")
  }
  
  output = list(
    variable       = variable_cut,
    bands          = bands,
    bands_count    = table(variable_cut),
    measure_value  = merging_summary[dim(merging_summary)[1],3],
    merging_raport = merging_summary
  )
  
  if( write_flag == TRUE & flag_optim == TRUE)
  {
    file_name = paste(var.name,"__pproc", sep = "")
    write.table( merging_summary, file = paste(write_file_path, file_name, ".csv", sep=""), sep = ";",dec = ",",row.names=F, na = "")
    #   file_name = paste(var.name,"__optim_partition", sep = "")
    #   write.table( bands, file = paste(write_file_path, file_name, ".csv", sep=""), sep = ";",dec = ",",row.names=F, na = "")
  }
  return( output )
  
}





#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
#---------------------------------------------------------------------- Historia zmian -------------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------#


## Nr wersji:          06
## Zakres zmian:
#   Korekta dla colnames(tab) gdy dim(tab) zwraca NULL


## Nr wersji:          05
## Zakres zmian:
#   Umo?liwienie obs?ugi targetu z przedzia?u [0,1] (uog?lnienie)
#############################################################################################################################################################################



#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
#---------------------------------------------------------- SPECYFIKACJA FUNKCJI 'bands_merging' ---------------------------------------------------------------------------#
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------#



# ----- CHARAKTERYSTYKA FUNKCJI --------------------------------------------------------------------------------------------------------------------------------------------#
# Funkcja scala kategorie skategoryzowanej, zadanej na wejsciu zmiennej "variable". Zmienna "variable" wraz z zmienn? "target" (z przedzia?u [0,1]) wyznaczaj? tabel? k x 2:
#
#         tagret=0 |target=1
#         - - - - - - - - -
# kat_1  |   n_01  |  n_11  |
#      - - - - - - - - - - -
# kat_2  |   n_02  |  n_12  |
#      - - - - - - - - - - -
#   .    |         |        |
#   .    |   ...   |  ...   |
#   .    |         |        |
#      - - - - - - - - - - -
# kat_k  |   n_0k  |  n_1k  |
#         - - - - - - - - -
# gdzie kat_1, ..., kat_k s? kategoriami zmiennej "variable". Jesli dla pewnych i, j zachodzi: n_ij <= threshold, to kategoria kat_i zostaje z??czona z kat_(i-1) lub z kat_(i+1).
# Nazwa nowej kategorii to po??czone znakiem "_" nazwy ??czonych kategorii.
# Mechanizm ??czenia kategorii:
# 1. Jako pierwsze ??czone s? kategorie 1, k. W dalszej kolejnosci: pierwsza z "wewn?trznych" kategorii.
# 2. Po??czenia (i-1,i) lub (i,i+1) s? klarowne, gdy scalane zostaj? z "s?siadem" kategorie 1 lub k. Po??czenia
#    (i-1,i) lub (i,i+1) w przypadku kategorii "wewn?trznych" wyjasnione s? w wierszach 74 -97.




# ----- OPIS ARGUMENT?W WEJSCIOWYCH ----------------------------------------------------------------------------------------------------------------------------------------#
# variable       -> zmienna, kt?rej kategorie podlegaj? scalaniu
# target         -> zmienna 0-1
# threshold      -> sta?a okreslaj?ca minimalny, dopuszczalny poziom licznosci kom?rki powyzszej tabeli
# labels_merged  -> data frame zawieraj?cy w kolejnych wierszach scalone juz kategorie ("variable" moze juz mie? niekt?re kategorie scalone z innych powod?w niz warunek
#                   "threshold")



# ----- OBIEKTY ZWRACANE ---------------------------------------------------------------------------------------------------------------------------------------------------#
#output - lista nast?puj?cych obiekt?w:
# variable       -> zmienna wejsciowa o scalonych kategoriach
# labels_merged  -> data frame zawieraj?cy kodowanie pozwalaj?ce uzyska? po??czenie kategorii zmiennej 'variable' kt?re ma miejsce w procesie ??czenia kategorii


#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------#



bands_merging = function(  variable, target, threshold = 0, labels_merged = data.frame() )
{
  variable        = factor(variable) #, exclude = NULL)
  tapp            = tapply(target, variable, FUN = "mean")
  ind_ordered     = order( as.numeric(tapp) )
  ####tab             = table(variable, target)[ind_ordered,]
  tab             = as.table( cbind(tapply(1-target, variable, FUN = "sum"), tapply(target, variable, FUN = "sum")) )[ind_ordered,] #uog?lnienie
  colnames(tab)   = c(0,1)
  tab[which(is.na(tab))] = 0 
  ind_under_thres = which(tab <= threshold, arr.ind = TRUE)
  i               = length(levels(variable))
  loop_end        = FALSE
  column          = 0
  row             = 0
  if( dim(labels_merged)[1] * dim(labels_merged)[2] == 0 )
  {
    labels_merged        = as.data.frame( matrix( NA, nrow = length(levels(variable)), ncol = length(levels(variable))+1 ) )
    labels_merged[,1]    = levels(variable)[ind_ordered]
    names(labels_merged) = c(paste("atr",1:length(levels(variable)),sep=""),"label")
    #   names(labels_merged) = levels(variable)[ind_ordered]
  }
  
  while( loop_end == FALSE )
  {
    if( length(ind_under_thres) > 0 )
    {
      if( is.element( 1, ind_under_thres[,1]) == TRUE )                      # czy pusty element pierwszego wiersza
      {
        row = 1
        row_merge = 2
      }
      else if( is.element( dim(tab)[1], ind_under_thres[,1]) == TRUE )  # czy pusty element ostatniego wiersza
      {
        row = dim(tab)[1]
        row_merge = dim(tab)[1] - 1
      }
      else                                                                 # pusty element wewn?trznego wiersza
      {
        row    = ind_under_thres[1,1]
        column = ind_under_thres[1,2]
        #       row_merge = row + c(-1,1)[sample(2,1)]
        s11 = sum(tab[c(row-1,row),1])
        s12 = sum(tab[c(row,row+1),1])
        s21 = sum(tab[c(row-1,row),2])
        s22 = sum(tab[c(row,row+1),2])
        if( s11 <= s12 & s21 <= s22 )
        {  row_merge = row - 1 }
        else if( s11 > s12 & s21 > s22 )
        {  row_merge = row + 1 }
        else if( sum(tab[,1]) > sum(tab[,2]) )
        {
          if( s21 <= s22 )
          { row_merge = row - 1 }
          else
          { row_merge = row + 1 }
        }
        else
        {
          if( s11 <= s12 )
          { row_merge = row - 1 }
          else
          { row_merge = row + 1 }
        }
      }
      
      if( row < row_merge )
      {
        ind_level_merge1 = which( levels(variable) == dimnames(tab)[[1]][row] )
        ind_level_merge2 = which( levels(variable) == dimnames(tab)[[1]][row_merge] )
      }
      else
      {
        ind_level_merge1 = which( levels(variable) == dimnames(tab)[[1]][row_merge] )
        ind_level_merge2 = which( levels(variable) == dimnames(tab)[[1]][row] )
      }
      ind1 = which( levels(variable)[ind_ordered] == levels(variable)[ind_level_merge1] )
      ind2 = which( levels(variable)[ind_ordered] == levels(variable)[ind_level_merge2] )
      if( ind1 < ind2 )
      {
        ind3                                    = which( is.na(labels_merged[ind1,]) )[1]
        ind4                                    = which( is.na(labels_merged[ind2,]) )[1]
        labels_merged[ ind1, ind3:(ind3+ind4-2)] = labels_merged[ ind2, 1:(ind4-1)]
        #  names(labels_merged)[ind1]              = paste( names(labels_merged)[ind1], names(labels_merged)[ind2], sep="_")
        #  labels_merged                           = as.data.frame(labels_merged[,-(ind2)])
        name_tmp                                = paste( labels_merged$label[ind1], labels_merged$label[ind2], sep="_")
        labels_merged                           = as.data.frame(labels_merged[-(ind2),])
        labels_merged$label[ind1]               = name_tmp
        
      }
      else
      {
        ind3                                    = which( is.na(labels_merged[ind2,]) )[1]
        ind4                                    = which( is.na(labels_merged[ind1,]) )[1]
        labels_merged[ind2, ind3:(ind3+ind4-2)] = labels_merged[ind1, 1:(ind4-1)]
        #       names(labels_merged)[ind2]              = paste( names(labels_merged)[ind2], names(labels_merged)[ind1], sep="_")
        #       labels_merged                           = labels_merged[,-(ind1)]
        name_tmp                                = paste( labels_merged$label[ind2], labels_merged$label[ind1], sep="_")
        labels_merged                           = as.data.frame(labels_merged[-(ind1),])
        labels_merged$label[ind2]               = name_tmp
      }
      
      label_new                          = paste( levels(variable)[ind_level_merge1], levels(variable)[ind_level_merge2], sep = "_")
      levels(variable)[ind_level_merge1] = label_new
      levels(variable)[ind_level_merge2] = label_new
      
      tapp                               = tapply(target, variable, FUN = "mean")
      ind_ordered                        = order( as.numeric(tapp) )
      ####tab                                = table(variable, target)[ind_ordered,]
      tab                                = as.table( cbind(tapply(1-target, variable, FUN = "sum"), tapply(target, variable, FUN = "sum")) )[ind_ordered,] #uog?lnienie
      if( !is.null(dim(tab)) )
      {colnames(tab)                     = c(0,1)
      }else
      {names(tab)                        = c(0,1)}
      tab[which(is.na(tab))]             = 0 
      ind_under_thres                    = which(tab <= threshold, arr.ind = TRUE)
    }
    else
    { loop_end = TRUE }
    i = i - 1
    if( i == 0 & loop_end == FALSE)
    { stop("Error during merging process") }
  }
  output = list(variable = variable, labels_merged = labels_merged)
  return(output)
}



#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
#---------------------------------------------------------------------- Historia zmian -------------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------#


## Nr wersji:          03
## Zakres zmian:
# Zwracanie etykiet do zakodowanych przedzia??w zmiennej zamiast WoE, gdy argument flag_label = TRUE (domy?lnie flag_label = FALSE aby zachowa? zgodno?? z poprzedni? wersj?)
#############################################################################################################################################################################





#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
#---------------------------------------------------------- SPECYFIKACJA FUNKCJI 'woe_coding' ------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------#


# ----- CHARAKTERYSTYKA FUNKCJI --------------------------------------------------------------------------------------------------------------------------------------------#
# Funkcja koduje zmienn? na wartosci WoE. Proces kodowania zalezy od typu zmiennej (ci?g?a/dyskretna).



# ----- ARGUMENTY WEJSCIOWE ------------------------------------------------------------------------------------------------------------------------------------------------#
# card      -> karta determinuj?ca kod WoE
# flag_cont -> wartos? logiczna TRUE/FALSE okreslaj?ca typ zmiennej (ci?g?a/dyskretna). ramka danych zawieraj?ca zmienne do zakodowania
# variable  -> zmienna kt?ra ma by? przekodowana na WoE
# dict      -> s?ownik mapuj?cy kategorie zmiennej dyskretnej(na inne kategorie przedstawione w s?owniku)



# ----- OBIEKTY ZWRACANE ---------------------------------------------------------------------------------------------------------------------------------------------------#
# var_woe -> zakodowana zmienna




#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------#




woe_coding = function(
  card                     ,
  flag_cont  = TRUE        ,
  variable                 ,
  dict       = data.frame(),
  flag_label = FALSE
)
{
  n           = dim(card)[1]
  k           = dim(card)[2]
  var_woe     = vector("numeric", length = length(variable))*NA
  var_label   = vector("character", length = length(variable))#*NA
  label_coded = vector("integer", length=n)
  flag_dict   = FALSE
  if( dim(dict)[1]*dim(dict)[1] > 0 )
  { flag_dict = TRUE  }
  
  par_return  = ifelse( flag_cont, 6, k)
  if( flag_label == TRUE )
  { par_return = par_return-1 }
  
  
  if( flag_cont == TRUE )
  {
    # -Inf and +Inf recognizing
    ind_inf_minus         = which( is.na( card[,2] ) * (card[,1]==0) == 1 )
    ind_inf_plus          = which( is.na( card[,3] ) * (card[,4]==0) == 1 )
    card[ind_inf_minus,2] = -Inf
    card[ind_inf_plus,3]  =  Inf
    
    len_variable  = length(variable)
    len_card      = length(card[,1]) # = dim(card)[1]
    
    lower_bounds  = matrix( rep(card[,2], each = len_variable), len_variable, len_card )
    upper_bounds  = matrix( rep(card[,3], each = len_variable), len_variable, len_card )
    variable_rep  = matrix( rep(variable, times = len_card), len_variable, len_card )
    
    lower_g       = lower_bounds < variable_rep
    lower_geq     = lower_bounds <= variable_rep
    upper_l       = upper_bounds > variable_rep
    upper_leq     = upper_bounds >= variable_rep
    
    ind_lower_geq = which( matrix( rep(card[,1], each = len_variable), len_variable, len_card )==1, arr.ind = TRUE )
    ind_upper_leq = which( matrix( rep(card[,4], each = len_variable), len_variable, len_card )==1, arr.ind = TRUE )
    
    lower         = lower_g
    upper         = upper_l
    
    if( is.matrix(ind_lower_geq) )
    {  lower[ind_lower_geq] = lower_geq[ind_lower_geq]  }
    
    if( is.matrix(ind_upper_leq) )
    { upper[ind_upper_leq] = upper_leq[ind_upper_leq]  }
    
    matrix_score_ind = lower*upper
    
    # BLANK (NA) recognizing
    ind_NoData = which( is.na( card[,2] ) * is.na( card[,3] ) * (card[,1]==1) * (card[,4]==1) == 1 )
    if( length(ind_NoData) > 0 )
    {
      ind_NA = which(is.na(variable))
      if( length(ind_NA) > 0 )
      {  matrix_score_ind[ind_NA, ind_NoData] = 1 }
      
    }
    ind_score              = which(matrix_score_ind == 1, arr.ind = TRUE)#[2]
    var_woe[ind_score[,1]] = card[ind_score[,2],par_return]
    #  ind_occur              = which(is.element(card[,5],unique(var_woe)))
    #  var_woe                = factor(var_woe)
    #  levels(var_woe)[1:length(ind_occur)] = card[ind_occur[order(card[ind_occur,6])],5]
  }
  #======================================================================================================================================#
  #============================= DISCRET VARIABLE =======================================================================================#
  #======================================================================================================================================#
  else
  {
    variable                        = factor(variable, exclude = NULL)
    ind_na                          = which(is.na(levels(variable)))
    levels(variable)[ind_na]        = "blank"
    if( flag_dict == TRUE )
    {
      dict$code                       = factor(dict$code, exclude = NULL)
      ind_na_dict                     = which( is.na(levels(dict$code)) )
      levels(dict$code)[ind_na_dict]  = "blank"
      dict_sorted                     = dict[order(dict$code),]
      l1                              = length(levels(dict$code))
      l2                              = length(levels(variable))
      matrix_ind                      = which( matrix( rep(levels(variable), times = l1 ), l2, l1 ) == matrix( rep(levels(dict$code), each = l2 ), l2, l1 ), arr.ind = TRUE)
      levels(variable)[matrix_ind[,1]]= dict_sorted$label[matrix_ind[,2]]
    }
    for( i in 1:n )
    {
      ind_aux            = ifelse( length( which(is.na(card[i,])) ) == 0, k-2, which(is.na(card[i,]))[1]-1 )
      ind_label          = which(is.element(variable,card[i,1:ind_aux]))
      var_woe[ind_label] = card[i,par_return]
      #   if( length(ind_label) > 0 )
      #   {  label_coded[i] = 1  }
      
    }
    # var_woe                                          = factor(var_woe, exclude = NULL)
    # levels(var_woe)[1:length(which(label_coded==1))] = c(card[order(card[which(label_coded==1),k]),k-1])
    # if( length(is.na(var_woe)) > 0)
    # { warning(paste(paste("Some label(s) haven't been met:"),paste(card[which(label_coded==0),k], sep=";")))  }
  }
  return(var_woe)
}
