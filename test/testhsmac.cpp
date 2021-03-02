//*****************************************************************************
//	Name		:	test/testhsmac.cpp
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

	//----------配列の確保----------
	double *p = new double[(nx + 2)*(ny + 2)];
    double *u = new double[(nx + 1)*(ny + 2)], *utilda = new double[(nx + 1)*(ny + 2)];
	double *v = new double[(nx + 2)*(ny + 1)], *vtilda = new double[(nx + 2)*(ny + 1)];
	
	//----------配列の初期化----------
	cout << "Initialized\n";
	for (int i = 0; i < (nx + 2)*(ny + 2); i++) {
		p[i] = 0.0;
	}
	for (int i = 0; i < (nx + 1)*(ny + 2); ++i) {
		u[i] = 0.0;
		utilda[i] = 0.0;
	}
	for (int i = 0; i < (nx + 2)*(ny + 1); ++i) {
		v[i] = 0.0;
		vtilda[i] = 0.0;
	}

	//----------HSMAC法の時間発展計算----------
	cout << "Culculate N-S Equation\n";
	for (int t = 0; t < tstep; t++) {
		cout << "\nt = " << (double)t*dt;
		//.....uの境界条件を設定.....
		for (int i = 1; i < nx; i++) {
			u[IndexU(i, 0)] = -u[IndexU(i, 1)];
			u[IndexU(i, ny + 1)] = 2.0 - u[IndexU(i, ny)];
		}

		//.....vの境界条件を設定.....
		for (int j = 1; j < ny; j++) {
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
				double udvdx = ((u[IndexU(i, j + 1)] + u[IndexU(i, j)])*(v[IndexV(i + 1, j)] + v[IndexV(i, j)]) - (u[IndexU(i - 1, j + 1)] + u[IndexU(i - 1, j)])*(v[IndexV(i,j)] + v[IndexV(i - 1, j)]))/(4.0*dx);
				double vdvdy = (pow(v[IndexV(i, j + 1)] + v[IndexV(i, j)], 2.0) - pow(v[IndexV(i, j)] + v[IndexV(i, j - 1)], 2.0))/(4.0*dy);
				double dpdy = (p[IndexP(i, j + 1)] - p[IndexP(i, j)])/dy;
				double d2vdx2 = (v[IndexV(i + 1, j)] - 2.0*v[IndexV(i, j)] + v[IndexV(i - 1, j)])/pow(dx, 2.0);
				double d2vdy2 = (v[IndexV(i, j + 1)] - 2.0*v[IndexV(i, j)] + v[IndexV(i, j - 1)])/pow(dy, 2.0);

				vtilda[IndexV(i, j)] = v[IndexV(i, j)] - dt*(udvdx + vdvdy + dpdy - (d2vdx2 + d2vdy2)/re);
			}
		}

		//.....D→0となるまで反復計算.....
		for (int k = 0; k < nk; k++) {
			double dmax = 0.0;
			//内部点の計算
			for (int i = 1; i < nx + 1; i++) {
				for (int j = 1; j < ny + 1; j++) {
					double larged = (utilda[IndexU(i, j)] - utilda[IndexU(i - 1, j)])/dx + (vtilda[IndexV(i, j)] - vtilda[IndexV(i, j - 1)])/dy;
					double deltap = -omega*larged/(2.0*dt*(pow(1.0/dx, 2.0) + pow(1.0/dy, 2.0)));
					if (fabs(deltap) > dmax) {
						dmax = fabs(deltap);
					}

					//.....pの更新.....
					p[IndexP(i, j)] += deltap;

					//.....utildaの更新.....
					if (i != 1) {
						utilda[IndexU(i - 1, j)] += -dt*deltap/dx;
					}
					if (i != nx) {
						utilda[IndexU(i, j)] += dt*deltap/dx;
					}

					//.....vtildaの更新.....
					if (j != 1) {
						vtilda[IndexV(i, j - 1)] += -dt*deltap/dy;
					}
					if (j != ny) {
						vtilda[IndexV(i, j)] += dt*deltap/dy;
					}
				}
			}			

			//収束判定
			if (dmax < eps) {
				cout << "\tConverged at k = " << k << "\tand dmax = " << dmax;
				break;
			}
		}	

		//.....uの値を更新.....
		for (int i = 1; i < nx; i++) {
			for (int j = 1; j < ny + 1; j++) {
				u[IndexU(i, j)] = utilda[IndexU(i, j)];
			}
		}

		//.....vの値を更新.....
		for (int i = 1; i < nx + 1; i++) {
			for (int j = 1; j < ny; j++) {
				v[IndexV(i, j)] = vtilda[IndexV(i, j)];
			}
		}
	}

	//----------結果の出力----------
	std::ofstream fout("result/hsmac.vtk");
    
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
    delete[] p, u, utilda, v, vtilda;

	return 0;
}