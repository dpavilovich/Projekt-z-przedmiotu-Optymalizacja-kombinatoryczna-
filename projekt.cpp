#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <ctime>
#include <vector>
#include <set>
#include <cstring>

using namespace std;


void loadDNAFromFile();
void generateInstanceManually();
void displayLastInstance();
void generateSpectrum(const string& dna, int k);
vector<string> generateReconstructedSpectrumWithDeltaK(const string& dna, int n, int k, int delta_k);
int totalLevenshteinDistance(const vector<string>& spectrum1, const vector<string>& spectrum2);
int calculateSimilarityPercentage(const vector<string>& original, const vector<string>& reconstructed);
int lastLevenshteinResult = 0;  
int lastSimilarityPercentage = 0; 




vector<string> reconstructedSpectrum;

vector<string> spectrum;




string lastGeneratedDNA = "";
int lastN = 0, lastK = 0, lastDeltaK = 0, lastLpoz = 0, lastLneg = 0;



void displayLastInstance() {
    if (lastGeneratedDNA.empty()) {
        cout << "Nie ma zapisanej instancji!" << endl;
        return;
    }

    cout << "\nOstatnia wygenerowana instancja:\n";
    cout << "DNA: " << lastGeneratedDNA << "\n";
    cout << "Parametry: n=" << lastN << ", k=" << lastK
        << ", delta_k=" << lastDeltaK << ", l_poz=" << lastLpoz
        << ", l_neg=" << lastLneg << "\n";
}


void loadDNAFromFile() {
    cout << "Wczytywanie DNA z pliku..." << endl;
    string filename;
    cout << "Podaj nazwę pliku: ";
    cin >> filename;

    ifstream file(filename);
    if (file.is_open()) {
        string dna;
        getline(file, dna);  
        file.close();
        cout << "Wczytane DNA: " << dna << endl;

        int n = dna.length();
        int k;
        int delta_k;

        
        cout << "Podaj długość k (długość oligonukleotydu): ";
        cin >> k;

        
        cout << "Podaj wartość delta_k (maksymalna zmiana długości oligonukleotydu): ";
        cin >> delta_k;

       
        generateSpectrum(dna, k);

        
        vector<string> reconstructedSpectrum = generateReconstructedSpectrumWithDeltaK(dna, n, k, delta_k);

        
        cout << "\nZrekonstruowane spektrum z delta_k:" << endl;
        for (const string& oligo : reconstructedSpectrum) {
            cout << oligo << endl;
        }

        
        char runNaive;
        cout << "\nCzy chcesz uruchomić algorytm naiwny? (T/N): ";
        cin >> runNaive;
        if (runNaive == 'T' || runNaive == 't') {
            int levenshteinResult = lastLevenshteinResult = totalLevenshteinDistance(spectrum, reconstructedSpectrum);
            cout << "Miara Levenshteina między spektrum oryginalnym a zrekonstruowanym: " << levenshteinResult << endl;
            
            int similarityPercentage = lastSimilarityPercentage = calculateSimilarityPercentage(spectrum, reconstructedSpectrum);

            
            cout << "Procent zgodności między spektrum oryginalnym a zrekonstruowanym: " << similarityPercentage << "%" << endl; 
        }
       
        char saveToFile;
        cout << "\nCzy chcesz zapisać wyniki do pliku? (T/N): ";
        cin >> saveToFile;
        if (saveToFile == 'T' || saveToFile == 't') {
            
            string fileName;
            cout << "Podaj nazwę pliku, do którego chcesz zapisać wyniki (np. wyniki.txt): ";
            cin >> fileName;

           
            ofstream outFile(fileName);
            if (outFile.is_open()) {
                
                outFile << "Generowane DNA: " << dna << endl;
                outFile << "\nSpektrum oryginalne:" << endl;
                for (const string& oligo : spectrum) {
                    outFile << oligo << endl;
                }

                outFile << "\nZrekonstruowane spektrum z delta_k:" << endl;
                for (const string& oligo : reconstructedSpectrum) {
                    outFile << oligo << endl;
                }

                outFile << "\nMiara Levenshteina między spektrum oryginalnym a zrekonstruowanym: " << lastLevenshteinResult << endl;
                outFile << "Procent zgodności między spektrum oryginalnym a zrekonstruowanym: " << lastSimilarityPercentage << "%" << endl;
                outFile << "\nk: " << k << ", delta_k: " << delta_k << endl;

                
                outFile.close();
                cout << "Wyniki zostały zapisane do pliku " << fileName << endl;
            }
            else {
                cout << "Nie udało się otworzyć pliku do zapisu." << endl;
            }
        }
    }
    else {
        cout << "Nie udało się otworzyć pliku!" << endl;
    }
}


char generateRandomNucleotide() {
    int r = rand() % 4;
    switch (r) {
    case 0: return 'A';
    case 1: return 'C';
    case 2: return 'G';
    case 3: return 'T';
    default: return 'A';
    }
}

string generateRandomDNA(int length) {
    string dna = "";
    for (int i = 0; i < length; i++) {
        dna += generateRandomNucleotide();  
    }
    return dna;
}

void generateSpectrum(const string& dna, int k) {
    int n = dna.length();

    
    vector<string> spectrum;
    for (int i = 0; i <= n - k; i++) {
        spectrum.push_back(dna.substr(i, k));  
    }

    
    cout << "Spektrum (oligonukleotydy o długości " << k << "):" << endl;
    for (const string& seq : spectrum) {
        cout << seq << endl;
    }
}

string generateRandomPositiveError(const string& oligo, int delta_k) {
    
    int k = oligo.length();
    int newLength = (rand() % 2 == 0) ? k + delta_k : k - delta_k;  

    
    newLength = max(1, min(newLength, (int)oligo.size()));  

    
    return oligo.substr(0, newLength);
}

vector<string> generatePositiveError(const string& oligo) {
    string error1 = oligo;
    string error2 = oligo;

    
    int lastIndex = error1.length() - 1;
    error1[lastIndex] = generateRandomNucleotide();

    
    int middleIndex = error1.length() / 2;
    if (error1.length() % 2 == 0) {
        middleIndex -= 1; 
    }
    error2[middleIndex] = generateRandomNucleotide();

    
    vector<string> positiveErrors = { error1, error2 };
    return positiveErrors;
}




vector<string> generatePositiveErrors(vector<string>& spectrum, int l_poz, int probablePositive, int delta_k) {
    vector<string> positiveErrors;

    
    if (l_poz > 0) {
        
        int numPositiveErrors = max(l_poz, 10);

        
        for (int i = 0; i < numPositiveErrors && i < spectrum.size(); i++) {
            string& oligo = spectrum[i]; 

            
            if (probablePositive == 0) {
                string positiveError = generateRandomPositiveError(oligo, delta_k);
                positiveErrors.push_back(positiveError);  
                oligo = positiveError;  
            }
            
            else if (probablePositive == 1) {
                vector<string> positiveError = generatePositiveError(oligo);
                positiveErrors.push_back(positiveError[0]); 
                positiveErrors.push_back(positiveError[1]); 
                
                oligo = positiveError[0];
            }
        }
    }

    return positiveErrors;  
}


vector<string> generateNegativeErrors(vector<string>& spectrum, int l_neg, int delta_k) {
    vector<string> negativeErrors;

    
    for (int i = 0; i < l_neg; i++) {
        int index = rand() % spectrum.size(); 
        string oligo = spectrum[index];

        
        int mutationIndex = rand() % oligo.length();
        char mutatedNucleotide = "ACGT"[rand() % 4];  
        oligo[mutationIndex] = mutatedNucleotide;

        
        if (delta_k > 0 && rand() % 2 == 0) {
            oligo.erase(mutationIndex, 1); 
        }

        negativeErrors.push_back(oligo);
        spectrum[index] = oligo;  
    }

    return negativeErrors;
}

vector<string> generateReconstructedSpectrumWithDeltaK(const string& dna, int n, int k, int delta_k) {
    vector<string> spectrum;
    srand(time(0)); 

    
    spectrum.push_back(dna.substr(0, k));

    
    for (int i = 1; i <= n - k; ++i) {
        
        int delta = rand() % (delta_k + 1);

        
        if (i >= n - k - 3) {
            
            spectrum.push_back(dna.substr(i, k));
        }
        else {
            if (delta == 0) {
               
                spectrum.push_back(dna.substr(i, k));
            }
            else {
                
                int sign = rand() % 2;  
                int newLength = k + (sign == 0 ? -delta : delta);

                
                if (i + newLength <= n && newLength > 0) {
                    spectrum.push_back(dna.substr(i, newLength));
                }
            }
        }
    }

    return spectrum;
}





void generateInstanceManually() {
    cout << "Generowanie ręczne..." << endl;



    char isStandardInstance;
    cout << "Czy instancja problemu ma być standardowa (T/N)? ";
    cin >> isStandardInstance;


   
    if (isStandardInstance == 'T' || isStandardInstance == 't') {
        int n = 400, k = 8, delta_k = 2, l_neg = 0, l_poz = 0, probablePositive = 0;
        bool repAllowed = true;

        
        cout << "Generowanie instancji z parametrami domyślnymi..." << endl;
        cout << "n: " << n << ", k: " << k << ", delta_k: " << delta_k << ", l_neg: " << l_neg
            << ", l_poz: " << l_poz << ", repAllowed: " << repAllowed << ", probablePositive: " << probablePositive << endl;

        
        string dna = lastGeneratedDNA = generateRandomDNA(n);
        cout << "\nWygenerowane DNA: " << dna << endl;
        generateSpectrum(dna, k);

        vector<string> spectrum;
        if (repAllowed) {
            for (int i = 0; i <= n - k; i++) {
                spectrum.push_back(dna.substr(i, k)); 
            }
        }
        else {
            set<string> uniqueOligos;
            for (int i = 0; i <= n - k; i++) {
                uniqueOligos.insert(dna.substr(i, k)); 
            }
            spectrum.assign(uniqueOligos.begin(), uniqueOligos.end());  
        }

        vector<string> positiveErrors = generatePositiveErrors(spectrum, l_poz, probablePositive, delta_k);

        
        cout << "Błędy pozytywne:" << endl;
        for (const string& error : positiveErrors) {
            cout << error << endl;
        }

        vector<string> negativeErrors = generateNegativeErrors(spectrum, l_neg, delta_k);

      
        cout << "Błędy negatywne:" << endl;
        for (const string& error : negativeErrors) {
            cout << error << endl;
        }

        
        cout << "\nZaktualizowane spektrum z błędami:" << endl;
        for (const string& oligo : spectrum) {
            cout << oligo << endl;
        }

        
        vector<string> reconstructedSpectrum = generateReconstructedSpectrumWithDeltaK(dna, n, k, delta_k);

        
        cout << "\nZrekonstruowane spektrum z delta_k:" << endl;
        for (const string& oligo : reconstructedSpectrum) {
            cout << oligo << endl;
        }

       
        char saveToFile;
        cout << "\nCzy chcesz zapisać instancję do pliku? (T/N): ";
        cin >> saveToFile;
        if (saveToFile == 'T' || saveToFile == 't') {
            
            ofstream outFile("instance.txt");
            if (outFile.is_open()) {
                outFile << "n: " << n << ", k: " << k << ", delta_k: " << delta_k << ", l_neg: " << l_neg
                    << ", l_poz: " << l_poz << ", repAllowed: " << repAllowed << ", probablePositive: " << probablePositive << endl;
                outFile << "\nWygenerowane DNA: " << dna << endl;
                outFile << "\nZaktualizowane spektrum:" << endl;
                for (const string& oligo : spectrum) {
                    outFile << oligo << endl;
                }
                outFile << "\nZrekonstruowane spektrum z delta_k:" << endl;
                for (const string& oligo : reconstructedSpectrum) {
                    cout << oligo << endl;
                }
                outFile.close();
                cout << "Instancja została zapisana do pliku instance.txt" << endl;
            }
            else {
                cout << "Nie udało się otworzyć pliku do zapisu." << endl;
            }
        }
       
        char runNaive;
        cout << "\nCzy chcesz uruchomić algorytm naiwny? (T/N): ";
        cin >> runNaive;
        if (runNaive == 'T' || runNaive == 't') {
            int levenshteinResult = lastLevenshteinResult = totalLevenshteinDistance(spectrum, reconstructedSpectrum);
            cout << "Miara Levenshteina między spektrum oryginalnym a zrekonstruowanym: " << levenshteinResult << endl;
            
            int similarityPercentage = lastSimilarityPercentage = calculateSimilarityPercentage(spectrum, reconstructedSpectrum);

            
            cout << "Procent zgodności między spektrum oryginalnym a zrekonstruowanym: " << similarityPercentage << "%" << endl;
        }

       
        char saveToFile2;
        cout << "\nCzy chcesz zapisać wyniki do pliku? (T/N): ";
        cin >> saveToFile2;
        if (saveToFile2 == 'T' || saveToFile2 == 't') {
            
            string fileName;
            cout << "Podaj nazwę pliku, do którego chcesz zapisać wyniki (np. wyniki.txt): ";
            cin >> fileName;

            
            ofstream outFile(fileName);
            if (outFile.is_open()) {
                
                outFile << "Generowane DNA: " << dna << endl;
                outFile << "\nSpektrum oryginalne:" << endl;
                for (const string& oligo : spectrum) {
                    outFile << oligo << endl;
                }

                outFile << "\nZrekonstruowane spektrum z delta_k:" << endl;
                for (const string& oligo : reconstructedSpectrum) {
                    outFile << oligo << endl;
                }

                outFile << "\nMiara Levenshteina między spektrum oryginalnym a zrekonstruowanym: " << lastLevenshteinResult << endl;
                outFile << "Procent zgodności między spektrum oryginalnym a zrekonstruowanym: " << lastSimilarityPercentage << "%" << endl;
                outFile << "\nk: " << k << ", delta_k: " << delta_k << endl;

                
                outFile.close();
                cout << "Wyniki zostały zapisane do pliku " << fileName << endl;
            }
            else {
                cout << "Nie udało się otworzyć pliku do zapisu." << endl;
            }
        }

        
        lastN = n;
        lastK = k;
        lastDeltaK = delta_k;
        lastLpoz = l_poz;
        lastLneg = l_neg;

        
        displayLastInstance();
    }
    else {
        cout << "Generowanie instancji ręczne..." << endl;
        
        int n = 400, k = 8, delta_k = 2, l_neg = 0, l_poz = 0, probablePositive = 0;
        bool repAllowed = true;

        
        cout << "Podaj wartosc n (dlugosc DNA, domyslnie 400): ";
        string input;
        cin.ignore();  
        getline(cin, input);  
        if (!input.empty()) {
            n = stoi(input);
            if (n < 300 || n > 700) {
                cout << "Wartość n musi być w zakresie 300-700. Używam domyślnej wartości 400." << endl;
                n = 400;
            }
        }

        cout << "Podaj wartość k (długość oligonukleotydów, domyślnie 8): ";
        getline(cin, input);
        if (!input.empty()) {
            k = stoi(input);
            if (k < 7 || k > 10) {
                cout << "Wartość k musi być w zakresie 6-10. Używam domyślnej wartości 8." << endl;
                k = 8;
            }
        }

        cout << "Podaj wartość delta_k (zakres zmian długości oligonukleotydów, domyślnie 2): ";
        getline(cin, input);
        if (!input.empty()) {
            delta_k = stoi(input);
            if (delta_k < 0 || delta_k > 2) {
                cout << "Wartość delta_k musi być w zakresie 0-2. Używam domyślnej wartości 2." << endl;
                delta_k = 2;
            }
        }

        
        cout << "Podaj wartość l_neg (liczba błędów negatywnych, domyślnie 0): ";
        getline(cin, input);
        if (!input.empty()) {
            l_neg = stoi(input);

            
            int max_neg_errors = n * 0.15;  
            if (l_neg < 10 || l_neg > max_neg_errors) {
                l_neg = 0; 
                cout << "Wartość l_neg jest niezgodna z warunkami. Ustawiam wartość domyślną l_neg = 0." << endl;
            }
        }

        cout << "Podaj wartość l_poz (liczba błędów pozytywnych, domyślnie 0): ";
        getline(cin, input);
        if (!input.empty()) l_poz = stoi(input);

        cout << "Czy dozwolone są powtórzenia w spektrum (1 = tak, 0 = nie, domyślnie 1): ";
        getline(cin, input);
        if (!input.empty()) repAllowed = (input == "1");

        cout << "Podaj wartość probablePositive (0 = losowe błędy, 1 = pozytywne błędy, domyślnie 0): ";
        getline(cin, input);
        if (!input.empty()) probablePositive = stoi(input);



        
        cout << "Generowanie instancji z parametrami: " << endl;
        cout << "n: " << n << ", k: " << k << ", delta_k: " << delta_k << ", l_neg: " << l_neg
            << ", l_poz: " << l_poz << ", repAllowed: " << repAllowed << ", probablePositive: " << probablePositive << endl;
        
        if (l_neg > 0) {
            repAllowed = true;  
        }

        
        string dna = lastGeneratedDNA = generateRandomDNA(n);
        cout << "\nWygenerowane DNA: " << dna << endl;
        generateSpectrum(dna, k);
        
        vector<string> spectrum;
        if (repAllowed) {
            for (int i = 0; i <= n - k; i++) {
                spectrum.push_back(dna.substr(i, k)); 
            }
        }
        else {
            set<string> uniqueOligos;
            for (int i = 0; i <= n - k; i++) {
                uniqueOligos.insert(dna.substr(i, k));  
            }
            spectrum.assign(uniqueOligos.begin(), uniqueOligos.end()); 
        }



        vector<string> positiveErrors = generatePositiveErrors(spectrum, l_poz, probablePositive, delta_k);

        
        cout << "Błędy pozytywne:" << endl;
        for (const string& error : positiveErrors) {
            cout << error << endl;
        }

        
        vector<string> negativeErrors = generateNegativeErrors(spectrum, l_neg, delta_k);

       
        cout << "Błędy negatywne:" << endl;
        for (const string& error : negativeErrors) {
            cout << error << endl;
        }

        
        cout << "\nZaktualizowane spektrum z błędami:" << endl;
        for (const string& oligo : spectrum) {
            cout << oligo << endl;
        }

        
        vector<string> reconstructedSpectrum = generateReconstructedSpectrumWithDeltaK(dna, n, k, delta_k);

        
        cout << "\nZrekonstruowane spektrum z delta_k:" << endl;
        for (const string& oligo : reconstructedSpectrum) {
            cout << oligo << endl;
        }



       
        cout << "\nCzy chcesz zapisać instancję do pliku? (1 = tak, 0 = nie): ";
        getline(cin, input);
        if (input == "1") {
            string fileName;
            cout << "Podaj nazwę pliku: ";
            getline(cin, fileName);

            ofstream outFile(fileName);  
            if (outFile.is_open()) {
               
                outFile << "Instancja DNA:\n";
                outFile << "n: " << n << ", k: " << k << ", delta_k: " << delta_k << ", l_neg: " << l_neg
                    << ", l_poz: " << l_poz << ", repAllowed: " << repAllowed << ", probablePositive: " << probablePositive << "\n";
                outFile << "Wygenerowane DNA: " << dna << "\n";
                outFile << "Spektrum DNA:\n";
                for (const string& oligo : spectrum) {
                    outFile << oligo << "\n";
                }
                outFile << "Błędy pozytywne:\n";
                for (const string& error : positiveErrors) {
                    outFile << error << "\n";
                }
                outFile << "Błędy negatywne:\n";
                for (const string& error : negativeErrors) {
                    outFile << error << "\n";
                }
                outFile << "\nZrekonstruowane spektrum z delta_k:";
                for (const string& oligo : reconstructedSpectrum) {
                    cout << oligo << endl;
                }

                outFile.close();
                cout << "Instancja została zapisana do pliku: " << fileName << endl;
            }
            else {
                cout << "Nie udało się otworzyć pliku do zapisu." << endl;
            }


        }
       
        char runNaive;
        cout << "\nCzy chcesz uruchomić algorytm naiwny? (T/N): ";
        cin >> runNaive;
        if (runNaive == 'T' || runNaive == 't') {
            
            int levenshteinResult = lastLevenshteinResult = totalLevenshteinDistance(spectrum, reconstructedSpectrum);
            cout << "Miara Levenshteina między spektrum oryginalnym a zrekonstruowanym: " << levenshteinResult << endl;
            
            int similarityPercentage = lastSimilarityPercentage = calculateSimilarityPercentage(spectrum, reconstructedSpectrum);

           
            cout << "Procent zgodności między spektrum oryginalnym a zrekonstruowanym: " << similarityPercentage << "%" << endl; 
        }

        
        char saveToFile;
        cout << "\nCzy chcesz zapisać wyniki do pliku? (T/N): ";
        cin >> saveToFile;
        if (saveToFile == 'T' || saveToFile == 't') {
           
            string fileName;
            cout << "Podaj nazwę pliku, do którego chcesz zapisać wyniki (np. wyniki.txt): ";
            cin >> fileName;

           
            ofstream outFile(fileName);
            if (outFile.is_open()) {
                
                outFile << "Generowane DNA: " << dna << endl;
                outFile << "\nSpektrum oryginalne:" << endl;
                for (const string& oligo : spectrum) {
                    outFile << oligo << endl;
                }

                outFile << "\nZrekonstruowane spektrum z delta_k:" << endl;
                for (const string& oligo : reconstructedSpectrum) {
                    outFile << oligo << endl;
                }

                outFile << "\nMiara Levenshteina między spektrum oryginalnym a zrekonstruowanym: " << lastLevenshteinResult << endl;
                outFile << "Procent zgodności między spektrum oryginalnym a zrekonstruowanym: " << lastSimilarityPercentage << "%" << endl;
                outFile << "\nk: " << k << ", delta_k: " << delta_k << endl;

               
                outFile.close();
                cout << "Wyniki zostały zapisane do pliku " << fileName << endl;
            }
            else {
                cout << "Nie udało się otworzyć pliku do zapisu." << endl;
            }
        }       
        lastN = n;
        lastK = k;
        lastDeltaK = delta_k;
        lastLpoz = l_poz;
        lastLneg = l_neg;


        
        displayLastInstance();
    }


}




int levenshteinDistance(const string& s1, const string& s2) {
    int len1 = s1.length();
    int len2 = s2.length();
    vector<vector<int>> dp(len1 + 1, vector<int>(len2 + 1));

    
    for (int i = 0; i <= len1; i++) {
        dp[i][0] = i;
    }
    for (int j = 0; j <= len2; j++) {
        dp[0][j] = j;
    }

    
    for (int i = 1; i <= len1; i++) {
        for (int j = 1; j <= len2; j++) {
            dp[i][j] = min({ dp[i - 1][j] + 1,           
                             dp[i][j - 1] + 1,           
                             dp[i - 1][j - 1] + (s1[i - 1] != s2[j - 1]) }); 
        }
    }

    return dp[len1][len2];
}


int totalLevenshteinDistance(const vector<string>& spectrum1, const vector<string>& spectrum2) {
    int totalDistance = 0;

    
    for (const string& oligo1 : spectrum1) {
        int minDist = INT_MAX;

        
        for (const string& oligo2 : spectrum2) {
            int dist = levenshteinDistance(oligo1, oligo2);
            minDist = min(minDist, dist);
        }

        
        totalDistance += minDist;
    }

    return totalDistance;
}

int calculateSimilarityPercentage(const vector<string>& original, const vector<string>& reconstructed) {
    int levenshteinDistanceResult = totalLevenshteinDistance(original, reconstructed);
    int maxLength = max(original.size(), reconstructed.size());  
    double similarityPercentage = (1 - static_cast<double>(levenshteinDistanceResult) / maxLength) * 100;
    return static_cast<int>(similarityPercentage);
}



int main() {
    int option;

    do {

        void displayLastInstance();
        
        cout << "\n--- MENU GLOWNE ---" << endl;
        cout << "1. Generator instancji" << endl;
        cout << "2. Algorytm naiwny" << endl;
        cout << "0. Zakoncz" << endl;
        cout << "Wybierz opcje: ";
        cin >> option;

        switch (option) {
        case 1:
            
            cout << "\n--- GENERATOR INSTANCJI ---" << endl;
            cout << "1. Wczytaj DNA z pliku" << endl;
            cout << "2. Generuj ręcznie" << endl;
            cout << "0. Powrót" << endl;
            cout << "Wybierz opcje: ";
            int subOption;
            cin >> subOption;

            switch (subOption) {
            case 1:
                loadDNAFromFile();  
                break;
            case 2:
                generateInstanceManually();  
                break;
            case 0:
                break;
            default:
                cout << "Niepoprawna opcja!" << endl;
                break;
            }
            break;

        case 2:
            
            cout << "Algorytm naiwny" << endl;
            cout << "Ostatnia miara Levenshteina: " << lastLevenshteinResult << endl;
            cout << "Ostatni procent zgodności: " << lastSimilarityPercentage << "%" << endl;
            break;

        case 0:
            cout << "Koniec programu." << endl;
            break;

        default:
            cout << "Niepoprawna opcja!" << endl;
            break;
        }
    } while (option != 0);

    return 0;
}
