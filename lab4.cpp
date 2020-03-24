#include <iostream>
#include <cmath>

using namespace std;

void mat(double **a, double *b, int n){ // Заполнение матрицы 
	for(int i = 0; i < n; i++){
		for(int j = 0; j < n; j++){
			if(i == j) a[i][j] = i+1;
			else if((i==(j+1)) || (i==(j-1))) a[i][j] = 1.;
			else a[i][j] = 0.;
		}
		b[i] = i+1; 
	}
}

void gauss(double **a, double *b, double *x, int m){ // Метод Гаусса
	for(int k = 1; k < m; k++){
		for(int j = k; j < m; j++){
			double h = a[j][k-1]/a[k-1][k-1];
			for(int i = 0; i < m; i++)
				a[j][i] -= h*a[k-1][i];
			b[j] -= h*b[k-1];
		}
		for(int i = m-1; i >= 0; i--){
			x[i] = b[i]/a[i][i];
			for(int c = m-1; c > i; c--)
				x[i] -= a[i][c]*x[c]/a[i][i];
		}
	}
}


void pvr(double **a, double *b, double *x, int m, double w, double eps){ // Метод ПВР
	double p[m];
	double delta;
	int k=0;
	do{
	    for(int i = 0; i < m; i++)
	        p[i] = x[i];

	    for(int i = 0; i < m; i++){
	        double var = 0.;
	        for(int j = 0; j < i; j++)
	            var += a[i][j]*x[j];

	        for(int j = i+1; j < m; j++)
	            var += a[i][j]*p[j];

	        x[i] = (1-w)*p[i] + w*(b[i] - var)/a[i][i];
	        
	    }
	    delta = 0;
	    for(int i = 0; i < m; i++)
	        delta += (x[i] - p[i])*(x[i] - p[i]);
	        k=k+1;
	} while(sqrt(delta) >= eps);
	cout<<"PVR iter="<<k<<endl;
}

void obr(double **a, int n){ // Поиск обратной матрицы
	double **b = new double*[n]; // Единичная матрица
	for(int i = 0; i < n; i++){
		b[i] = new double[n];
		b[i][i] = 1.;
	}
	double **ab = new double*[n];
	for(int i = 0; i < n; i++) ab[i] = new double[2*n];

	for(int i = 0; i < n; i++){
	    for(int j = 0; j < n; j++){
	        ab[i][j] = a[i][j];
	        ab[i][j + n] = b[i][j];
	    }
	}
	// Зануление нижнего левого угла
	for(int k = 0; k < n; k++){
	    for(int i = 0; i < 2*n; i++)
	        ab[k][i] = ab[k][i]/a[k][k];

		for(int i = k + 1; i < n; i++){
	        double K = ab[i][k]/ab[k][k];
			for(int j = 0; j < 2*n; j++)
	            ab[i][j] -= ab[k][j]*K;
	    }
	    for(int i = 0; i < n; i++)
	        for(int j = 0; j < n; j++)
	            a[i][j] = ab[i][j];
	}
	// Зануление верхнего правого угла
	for(int k = n - 1; k > -1; k--){
	    for(int i = 2*n - 1; i > -1; i--)
	        ab[k][i] = ab[k][i]/a[k][k];

	    for(int i = k - 1; i > -1; i--){
	        double K = ab[i][k]/ab[k][k];
	        for(int j = 2*n - 1; j > -1; j--)
	            ab[i][j] -= ab[k][j]*K;
	    }
	}

	for(int i = 0; i < n; i++)
	    for(int j = 0; j < n; j++)
	        a[i][j] = ab[i][j + n];
	delete[] b;
	delete[] ab;
}

double norm(double *a, int n){ // Функция нахождения нормы вектора
	double norm = 0.;
	for(int i = 0; i < n; i++)
		norm += a[i]*a[i];
	return sqrt(norm);
}

double eig(double **a, int n, double eps){ // Функция поиска максимального
	double x[n] = {};					   // собственного значения
	x[0] = 1.;						 	   // степенным методом 
	double y[n];
	double lc = 0.;
	double ln = 0.;
	do{
		lc = ln;
		for(int i = 0; i < n; i++){
			y[i] = 0.;
			for(int j = 0; j < n; j++)
				y[i] += a[i][j]*x[j];
		}
		ln = 0.;
		for(int i = 0; i < n; i++)
			ln += x[i]*y[i]/pow(norm(x, n),2);

		for(int i = 0; i < n; i++)
			x[i] = y[i]/norm(y, n);

	} while(abs(ln - lc) >= eps);
	return ln;
}

void pvr_dop(double *x, int m, double w, double eps){ // ПВР для больших матриц
	double p[m];								  	 
	double cov;
	int k=0;
	do{
	    for(int i = 0; i < m; i++)
	        p[i] = x[i];

	    for(int i = 0; i < m; i++){
	        double var = 0.;
	        for(int j = 0; j < i; j++){
				if(i == j) var += 10.*x[j];
				else if(i - 5 < j && j  < i + 5) var += x[j];
			}
	        for(int j = i+1; j < m; j++){
				if(i == j) var += 10.*p[j];
				else if(i - 5 < j && j  < i + 5) var += p[j];
			}
	        x[i] = (1-w)*p[i] + w*(i+1-var)/10.;
	    }
	    cov = 0.;
	    for(int i = 0; i < m; i++)
	        cov += (x[i] - p[i])*(x[i] - p[i]);
	        k++;
	} while(sqrt(cov) >= eps);
	cout<<"pvr na dop iter="<<k<<endl;
}

int main(){
	int n = 100;
	double **a = new double*[n];
	for(int i = 0; i < n; i++) a[i] = new double[n];
	double b[n];
	double x[n] = {};

	mat(a, b, n);

	pvr(a, b, x, n, 0.5, 1e-6);

	cout << "RESH PVR\n";
	for(int i = 0; i < n; i++)
		cout << x[i] << " ";
	cout << endl;
	
	// Вычисление вектора невязки
	double r[n] = {};
	for(int i = 0; i < n; i++){
		double s = 0.;
		for(int j = 0; j < n; j++)
			s += a[i][j]*x[j];
		r[i] = b[i] - s;
	}
	cout << "NORMA VECTORA NEVIAZKI  " << norm(r, n) << endl;

	gauss(a, b, x, n);

	cout << "\nRESHENIE GAUS\n";
	for(int i = 0; i < n; i++)
		cout << x[i] << " ";
	cout << endl;

	mat(a, b, n);
	// Вычисление вектора невязки
	for(int i = 0; i < n; i++){
		double s = 0.;
		for(int j = 0; j < n; j++)
			s += a[i][j]*x[j];
		r[i]= b[i] - s;
	}
	cout << "NORMA VECTORA NEVIAZKI " << norm(r, n) << endl;

	double mu = 1.; // Число обусловленности
	double ax_norm = 0.;
	for(int i = 0; i < n; i++){
		for(int j = 0; j < n; j++)
			ax_norm += pow(a[i][j]*x[j], 2);
	}
	mu *= sqrt(ax_norm)/pow(norm(x, n), 2);

	double inv_ax_norm = 0.;
	obr(a, n);
	for(int i = 0; i < n; i++)
		for(int j = 0; j < n; j++)
			inv_ax_norm += pow(a[i][j]*x[j], 2);
	mu *= sqrt(inv_ax_norm);

	cout << "\nCHISLO OBUSLOVLENNOSTI " << mu << endl;

	mat(a, b, n);

	// Дополнительные задания
	double l = eig(a, n, 1e-8);
	cout << "\nMAX SOBSTVENNOE ZNACH " << l << endl;
	obr(a, n);
	l = eig(a, n, 1e-10);
	cout << "MIN SOBSTVENNOE ZNACH " << 1./l << endl;


	cout << "\n RESHENIE MATR NE POMECHAETTSY\nVVEDITE RAZMER MATR n: ";
	cin >> n;
	double xn[n] = {};

	pvr_dop(xn, n, 1.5, 1e-4);

	double rn[n] = {};
	for(int i = 0; i < n; i++){
		double s = 0.;
		for(int j = 0; j < n; j++){
			if(i == j) s += 10.*xn[j];
			else if(i - 5 < j && j < i + 5) s+= xn[j];
		}
		rn[i] = i+1-s;
	}
	cout << "\nNORMA VEKTORA NESBAZKI " << norm(rn, n) << endl;

	delete[] a;
	return 0;
}
