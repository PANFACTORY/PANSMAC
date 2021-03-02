//*****************************************************************************
//	Name		:	test/testsmac.cpp
//	Author		:	Tanabe Yuta
//	Date		:	2021/03/02
//	Copyright	:	(C)2021 TanabeYuta
//*****************************************************************************

#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <fstream>

using namespace std;

int main() {
    int nx = 100, ny = 100, nk = 1000000;
    double tmax = 100.0, dt = 0.01, dx = 1.0/(double)nx, dy = 1.0/(double)ny, re = 1000.0, eps = 1e-8, omega = 2.0/(1.0 + sin(M_PI/(double)nx));
    int tstep = (int)(tmax/dt);

    //----------インデックス取得関数----------
    auto IndexP = [&](int _i, int _j) {
        return (ny + 2)*_i + _j;
    };

    auto IndexU = [&](int _i, int _j) {
        return (ny + 2)*_i + _j;
    };

    auto IndexV = [&](int _i, int _j) {
        return (ny + 1)*_i + _j;
    };

	//----------状態量を保存する配列----------
	double *p = new double[(nx + 2)*(ny + 2)], *phi = new double[(nx + 2)*(ny + 2)];
	double *u = new double[(nx + 1)*(ny + 2)], *utilda = new double[(nx + 1)*(ny + 2)];
	double *v = new double[(nx + 2)*(ny + 1)], *vtilda = new double[(nx + 2)*(ny + 1)];
	
	//----------配列の初期化----------
	cout << "Initialized\n";
	for (int i = 0; i < (nx + 2)*(ny + 2); ++i) {
		p[i] = 0.0;
		phi[i] = 0.0;
	}
	for (int i = 0; i < (nx + 1)*(ny + 2); ++i) {
		u[i] = 0.0;
		utilda[i] = 0.0;
	}
	for (int i = 0; i < (nx + 2)*(ny + 1); ++i) {
		v[i] = 0.0;
		vtilda[i] = 0.0;
	}

	//----------SMAC法の時間発展計算----------
	cout << "Culculate N-S Equation\n";
	for (int t = 0; t < tstep; t++) {
		cout << "\nt = " << (double)t*dt;

		//.....uの境界条件を設定.....
		for (int i = 0; i < nx + 1; i++) {
			u[IndexU(i, 0)] = -u[IndexU(i, 1)];
			u[IndexU(i, ny + 1)] = 2.0 - u[IndexU(i, ny)];
		}
 
		//.....vの境界条件を設定.....
		for (int j = 0; j < ny + 1; j++) {
			v[IndexV(0, j)] = -v[IndexV(1, j)];
			v[IndexV(nx + 1, j)] = -v[IndexV(nx, j)];
		}
 
		//.....utildaの計算.....
		for (int i = 1; i < nx; i++) {
			for (int j = 1; j < ny + 1; j++) {
				double ududx = (pow(u[IndexU(i + 1, j)] + u[IndexU(i, j)], 2.0) - pow(u[IndexU(i, j)] + u[IndexU(i - 1, j)], 2.0))/(4.0*dx);
				double vdudy = ((u[IndexU(i, j + 1)] + u[IndexU(i, j)])*(v[IndexV(i + 1, j)] + v[IndexV(i, j)]) - (u[IndexU(i, j)] + u[IndexU(i, j - 1)])*(v[IndexV(i + 1, j - 1)] + v[IndexV(i, j - 1)]))/(4.0*dy);
				double dpdx = (p[IndexP(i + 1, j)] - p[IndexP(i, j)])/dx;
				double d2udx2 = (u[IndexU(i + 1, j)] - 2.0*u[IndexU(i, j)] + u[IndexU(i - 1, j)])/pow(dx, 2.0);
				double d2udy2 = (u[IndexU(i, j + 1)] - 2.0*u[IndexU(i, j)] + u[IndexU(i, j - 1)])/pow(dy, 2.0);

				utilda[IndexU(i, j)] = u[IndexU(i, j)] - dt*(ududx + vdudy + dpdx - (d2udx2 + d2udy2)/re);
			}
		}

		//.....vtildaの計算.....
		for (int i = 1; i < nx + 1; i++) {
			for (int j = 1; j < ny; j++) {
				double udvdx = ((u[IndexU(i, j + 1)] + u[IndexU(i, j)])*(v[IndexV(i + 1, j)] + v[IndexV(i, j)]) - (u[IndexU(i - 1, j + 1)] + u[IndexU(i - 1, j)])*(v[IndexV(i, j)] + v[IndexV(i - 1, j)]))/(4.0*dx);
				double vdvdy = (pow(v[IndexV(i, j + 1)] + v[IndexV(i, j)], 2.0) - pow(v[IndexV(i, j)] + v[IndexV(i, j - 1)], 2.0))/(4.0*dy);
				double dpdy = (p[IndexP(i, j + 1)] - p[IndexP(i, j)])/dy;
				double d2vdx2 = (v[IndexV(i + 1, j)] - 2.0*v[IndexV(i, j)] + v[IndexV(i - 1, j)])/pow(dx, 2.0);
				double d2vdy2 = (v[IndexV(i, j + 1)] - 2.0*v[IndexV(i, j)] + v[IndexV(i, j - 1)])/pow(dy, 2.0);

				vtilda[IndexV(i, j)] = v[IndexV(i, j)] - dt*(udvdx + vdvdy + dpdy - (d2vdx2 + d2vdy2)/re);
			}
		}

		//.....φが収束するまで反復計算.....
		for (int k = 0; k < nk; k++) {
			double phimax = 0.0;
			//内部点の計算
			for (int i = 1; i < nx + 1; i++) {
				for (int j = 1; j < ny + 1; j++) {
					if (i == 1 && j == 1) {
						phi[IndexP(i, j)] = 0.0;
					} else {
						double tmpphi = phi[IndexP(i, j)];
						double larged = (utilda[IndexU(i, j)] - utilda[IndexU(i - 1, j)])/dx + (vtilda[IndexV(i, j)] - vtilda[IndexV(i, j - 1)])/dy;
						phi[IndexP(i, j)] = (1.0 - omega)*phi[IndexP(i, j)] + omega*pow(dx*dy, 2.0)/(2.0*(pow(dx, 2.0) + pow(dy, 2.0)))*((phi[IndexP(i + 1, j)] + phi[IndexP(i - 1, j)])/pow(dx, 2.0) + (phi[IndexP(i, j + 1)] + phi[IndexP(i, j - 1)])/pow(dy, 2.0) + larged);
						if (abs(phi[IndexP(i, j)] - tmpphi) > phimax) {
							phimax = abs(phi[IndexP(i, j)] - tmpphi);
						}
					}
				}
			}
			//仮想セルの計算
			for (int i = 1; i < nx + 1; i++) {
				phi[IndexP(i, 0)] = phi[IndexP(i, 1)];
				phi[IndexP(i, ny + 1)] = phi[IndexP(i, ny)];
			}
			for (int j = 1; j < ny + 1; j++) {
				phi[IndexP(0, j)] = phi[IndexP(1, j)];
				phi[IndexP(nx + 1, j)] = phi[IndexP(nx, j)];
			}
			//収束判定
			if (phimax < eps) {
				cout << "\tConverged at k = " << k;
				break;
			}
		}

		//.....uの値を更新.....
		for (int i = 1; i < nx; i++) {
			for (int j = 1; j < ny + 1; j++) {
				u[IndexU(i, j)] = utilda[IndexU(i, j)] + (phi[IndexP(i + 1, j)] - phi[IndexP(i, j)])/dx;
			}
		}

		//.....vの値を更新.....
		for (int i = 1; i < nx + 1; i++) {
			for (int j = 1; j < ny; j++) {
				v[IndexV(i, j)] = vtilda[IndexV(i, j)] + (phi[IndexP(i, j + 1)] - phi[IndexP(i, j)])/dy;
			}
		}

		//.....pの値の更新.....
		for (int i = 0; i < nx + 2; i++) {
			for (int j = 0; j < ny + 2; j++) {
				p[IndexP(i, j)] -= phi[IndexP(i, j)]/dt;
			}
		}
	}

	//----------結果の出力----------
	std::ofstream fout("result/smac.vtk");
    
    //  Add vtk header  
    fout << "# vtk DataFile Version 3.0" << std::endl;
    fout << "cavity" << std::endl;
    fout << "ASCII" << std::endl;
    fout << "DATASET\tSTRUCTURED_GRID" << std::endl;
    fout << "DIMENSIONS\t" << nx + 1 << "\t" << ny + 1 << "\t" << 1 << std::endl;

    //  Add point coordinates
    fout << "POINTS\t" << (nx + 1)*(ny + 1) << "\t" << "float" << std::endl;
    for (int j = 0; j < ny + 1; j++) {
        for (int i = 0; i < nx + 1; i++) {
            fout << i*dx << "\t" << j*dy << "\t" << 0.0 << std::endl;
        }
    }

    //  Export velocity
    fout << "CELL_DATA\t" << nx*ny << std::endl;
    fout << "VECTORS\tu\tfloat" << std::endl;
    for (int j = 0; j < ny; j++) {
        for (int i = 0; i < nx; i++) {
            fout << 0.5*(u[IndexU(i, j + 1)] + u[IndexU(i + 1, j + 1)]) << "\t" << 0.5*(v[IndexV(i + 1, j + 1)] + v[IndexV(i + 1, j)]) << "\t" << 0.0 << std::endl;
        }
    }

    //  Export pressure
    fout << "SCALARS\tp\tfloat" << std::endl;
    fout << "LOOKUP_TABLE\tdefault" << std::endl;
    for (int j = 1; j < ny + 1; j++) {
        for (int i = 1; i < nx + 1; i++) {
            fout << p[IndexP(i, j)] << std::endl;
        }
    }

	//----------メモリの開放----------
    delete[] p, phi, u, utilda, v, vtilda;
	
    return 0;
}