//Preston Peck
//CS 365
//September 14, 2017
//HW 2

#include <iostream>
#include <cmath>
using namespace std;

double price_from_yield(double F, double c, double y, int n, double B);
void price_from_yield(double F, double c, double y, int n, double &B, double &D_mac, double &D_mod);
int yield_from_price_bisection(double F, double c, int n, double B_market, double tol, int max_iter, double &y, int &num_iter);
int yield_from_price_NR(double F, double c, int n, double B_market, double tol, int max_iter, double &y, int &num_iter);
int yield_from_price_secant(double F, double c, int n, double B_market, double tol, int max_iter, double &y, int &num_iter);
int yield_from_price_fixpt(double F, double c, int n, double B_market, double tol, int max_iter, double &y, int &num_iter);

int main() {
    double F, c, B_market, tol, y = 0.0;
    int n, max_iter, num_iter = 0;
    cout << "ENTER VALUES" <<endl;
    cout << "Face Value: $";
    cin >> F;
    cout << "Coupon Payment: $";
    cin >> c;
    cout << "Payments: ";
    cin >> n;
    cout << "Market Value: $";
    cin >> B_market;
    cout << "Tolerance: ";
    cin >> tol;
    cout << "Max Iteration: ";
    cin >> max_iter;
    cout << endl;

    yield_from_price_bisection(F, c, n, B_market, tol, max_iter, y, num_iter);
    yield_from_price_NR(F, c, n, B_market, tol, max_iter, y, num_iter);
    yield_from_price_secant(F, c, n, B_market, tol, max_iter, y, num_iter);
    yield_from_price_fixpt(F, c, n, B_market, tol, max_iter, y, num_iter);
}

void price_from_yield(double F, double c, double y, int n, double &B, double &D_mac, double &D_mod) {
	double y_decimal = .01 * y;
	B = 0.0;
	D_mac = 0.0;
	D_mod = 0.0;
	
	for (int i = 1; i <= n; i++) {
		if (i < n) {
			B += (.5 * c) / pow(1 + (.5 * y_decimal), i);
		}
		
		else {
			B += (F + (.5 * c)) / pow(1 + (.5 * y_decimal), n);
		}
	}
	
	cout << "Maturity Value: $" << B << endl;
	
	for (int i = 1; i <= n; i++) {
		if (i < n) {
			D_mac += (.5 * i) * ((.5 * c) / pow(1 + (.5 * y_decimal), i));
		}
		
		else {
			D_mac += (.5 * i) * ((F + (.5 * c)) / pow(1 + (.5 * y_decimal), i));
			D_mac *= 1 / B;
			cout << "Macaulay Duration: " << D_mac << endl;
		}
	}
	
	D_mod = D_mac / (1 + (.5 * y_decimal));
	cout << "Modified Duration: " << D_mod << endl;
}

int yield_from_price_bisection(double F, double c, int n, double B_market, double tol, int max_iter, double &y, int &num_iter) {
    cout << "BISECTION" << endl;
    num_iter = 0;
    cout << "ITERATION " << 1 << endl;
    cout << "Low ";
    double y_low = 0.0;
    double By_low = price_from_yield(F, c, y_low, n, 0.0);

    if (abs(By_low - B_market) <= tol) {
    	y = y_low;
    	cout << "Yield Rate: " << y << '%' << '\n' << endl;
    	return 0;
    }

    if (By_low < B_market) {
    	y = 0;
    	cout << "Yield Rate: " << y << '%' << '\n' << endl;
    	return 1;
    }

    double y_high = 100.0;
    cout << "High ";
    double By_high = price_from_yield(F, c, y_high, n, 0.0);

    if (abs(By_high - B_market) <= tol) {
    	y = y_high;
    	cout << "Yield Rate: " << y << '%' << endl;
    	cout << endl << "Iterations: " << (num_iter + 1) << endl << endl;
    	return 0;
    }

    if (By_high > B_market) {
    	y = 0;
    	cout << "Yield Rate: " << y << '%' << endl;
    	cout << endl << "Iterations: " << (num_iter + 1) << '\n' << endl << endl;
    	return 1;
    }

    for (int i = 0; i < max_iter; ++i) {
        if (i >= 1) {
            cout << "ITERATION " << (i + 1) << endl;
        }

        y = (y_low + y_high) / 2.0;
        double B = price_from_yield(F, c, y, n, 0.0);
    
    	if (abs(B - B_market) <= tol || y_high - y_low <= tol) {
    		num_iter = i;
    		cout << "Yield Rate: " << y << '%' << endl;
    		cout << endl << "Iterations: " << (num_iter + 1) << endl << endl;
    		return 0;
    	}
    
    	else if (B > B_market) {
            cout << "Lower Limit: " << y << '%' << endl;
    		y_low = y;
    	}
    
    	else {
            cout << "Upper Limit: " << y << '%' << endl;
    		y_high = y;
    	}
    }

    num_iter = max_iter;
    y = 0;
    cout << "Yield Rate: " << y << '%' << endl;
    cout << endl << "Iterations: " << (num_iter + 1) << endl << endl;
    return 1;
}

int yield_from_price_NR(double F, double c, int n, double B_market, double tol, int max_iter, double &y, int &num_iter) {
    cout << "NEWTON RAPHSON" << endl;
    double f_iter[max_iter + 1];
	double y_iter[max_iter + 1];
	
	y_iter[0] = c;
	double B = 0.0;
	double D_mac = 0.0;
	double D_mod = 0.0;
	
	for (int i = 0; i < max_iter; ++i) {
        cout << "ITERATION " << (i + 1) << endl;
		price_from_yield(F, c, y_iter[i], n, B, D_mac, D_mod);
		f_iter[i] = B - B_market;
		
		if (abs(f_iter[i]) <= tol) {
			y = y_iter[i];
			cout << "Yield Rate: " << y << '%' << endl;
			num_iter = i;
			cout << endl << "Iterations: " << (num_iter + 1) << endl << endl;
			return 0;
		}
		
		double fprime = -0.01 * B * D_mod;
		
		if (fprime == 0.0) {
			y = 0.0;
			cout << "Yield Rate: " << y << '%' << endl;
			num_iter = 0;
			cout << endl << "Iterations: " << (num_iter + 1) << endl << endl;
			return 1;
		}
		
		double delta_y = f_iter[i] / fprime;
		
		if (abs(delta_y) <= tol) {
			y = y_iter[i];
			cout << "Yield Rate: " << y << '%' << endl;
			num_iter = i;
			cout << endl << "Iterations: " << (num_iter + 1) << endl << endl;
			return 0;
		}
	
		y_iter[i+1] = y_iter[i] - delta_y;
		cout << "Yield Rate: " << y_iter[i] << '%' << endl;
	}
	
	y = 0;
	cout << "Yield Rate: " << y << '%' <<endl;
	num_iter = max_iter;
	cout << endl << "Iterations: " << (num_iter + 1) << endl << endl;
	return 1;
}

int yield_from_price_secant(double F, double c, int n, double B_market, double tol, int max_iter, double &y, int &num_iter) {
    cout << "SECANT" << endl;
    double f_iter[max_iter + 1];
	double y_iter[max_iter + 1]; 
	
	y_iter[0] = c;
	y_iter[1] = c + 0.001;
	double B = 0.0;
	
	cout << "ITERATION " << 1 << endl;
	B = price_from_yield(F, c, y_iter[0], n, B);
	cout << "Yield Rate: " << y_iter[0] << '%' << endl;
    f_iter[0] = B - B_market;
    
    for (int i = 1; i < max_iter; ++i) {
    	cout << "ITERATION " << (i + 1) << endl;
	    B = price_from_yield(F, c, y_iter[i], n, B);
	    f_iter[i] = B - B_market;
    	double fprime = (f_iter[i] - f_iter[i-1]) / (y_iter[i] - y_iter[i-1]);
    	
    	if (fprime == 0.0) {
			y = 0.0;
			cout << "Yield Rate: " << y << '%' << endl;
			num_iter = 0;
			cout << endl << "Iterations: " << (num_iter + 1) << endl << endl;
			return 1;
		}
		
		double delta_y = f_iter[i] / fprime;
		
		if (abs(delta_y) <= tol) {
			y = y_iter[i];
			cout << "Yield Rate: " << y << '%' << endl;
			num_iter = i;
			cout << endl << "Iterations: " << (num_iter + 1) << endl << endl;
			return 0;
		}
	
		y_iter[i+1] = y_iter[i] - delta_y;
		cout << "Yield Rate: " << y_iter[i+1] << '%' << endl;
	}
	
	y = 0;
	cout << "Yield Rate: " << y << '%' << endl;
	num_iter = max_iter;
	cout << endl << "Iterations: " << (num_iter + 1) << endl << endl;
	return 1;
}

int yield_from_price_fixpt(double F, double c, int n, double B_market, double tol, int max_iter, double &y, int &num_iter) {
    cout << "FIXED POINT" << endl;
    double y_iter[max_iter + 1];
	y_iter[0] = c;
	double scale = 10.0;
	double B = 0.0;
	
	for (int i = 0; i < max_iter; ++i) {
		cout << "ITERATION " << (i + 1) << endl;
		B = price_from_yield(F, c, y_iter[i], n, B);
        y_iter[i+1] = y_iter[i] + scale * (B - B_market) / B_market;
        
        if (abs(B - B_market) <= tol || abs(y_iter[i+1] - y_iter[i]) <= tol) {
        	y = y_iter[i];
        	cout << "Yield Rate: " << y << '%' << endl;
        	num_iter = i;
        	cout << endl << "Iterations: " << (num_iter + 1) << endl << endl;
        	return 0;
        }
        
        cout << "Yield Rate: " << y_iter[i] << '%' << endl;
	}
	
	y = 0;
	cout << "Yield Rate: " << y << '%' << endl;
	num_iter = max_iter;
	cout << endl << "Iterations: " << (num_iter + 1) << endl << endl;
	return 1;
}

double price_from_yield(double F, double c, double y, int n, double B) {
	double y_decimal = .01 * y;
	B = 0.0;
	
	for (int i = 1; i <= n; i++) {
		if (i < n) {
			B += (.5 * c) / pow(1 + (.5 * y_decimal), i);
		}
		
		else {
			B += (F + (.5 * c)) / pow(1 + (.5 * y_decimal), n);
		}
	}
	
	cout << "Maturity Value: $" << B << endl;
	return B;
}