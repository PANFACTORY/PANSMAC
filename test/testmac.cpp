//*****************************************************************************
//	Name		:	test/testmac.cpp
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
    double tmax = 100.0, dt = 0.001, dx = 1.0/(double)nx, dy = 1.0/(double)ny, re = 1000.0, eps = 1e-8, omega = 2.0/(1.0 + sin(M_PI/(double)nx));
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
    double *u = new double[(nx + 1)*(ny + 2)], *up1 = new double[(nx + 1)*(ny + 2)];
    double *v = new double[(nx + 2)*(ny + 1)], *vp1 = new double[(nx + 2)*(ny + 1)];
	
	//----------配列の初期化----------
	cout << "Initialized\n";
	for (int i = 0; i < (nx + 2)*(ny + 2); i++) {
		p[i] = 0.0;
	}
	for (int i = 0; i < (nx + 1)*(ny + 2); i++) {
        u[i] = 0.0;
        up1[i] = 0.0;
	}
	for (int i = 0; i < (nx + 2)*(ny + 1); i++) {
		v[i] = 0.0;
		vp1[i] = 0.0;
	}
	
	//----------N-S方程式を時間発展的に解く----------
	cout << "Culculate N-S Equation\n";
	for (int t = 1; t <= tstep; t++) {
		cout << "\nt = " << (double)t*dt;

		//.....速度の境界条件を設定.....
		//上下壁仮想点
		for (int i = 0; i < nx + 1; i++) {
			u[IndexU(i, 0)] = -u[IndexU(i, 1)];
			u[IndexU(i, ny + 1)] = 2.0 - u[IndexU(i, ny)];
		}
		//左右壁仮想点
		for (int j = 0; j < ny + 1; j++) {
			v[IndexV(0, j)] = -v[IndexV(1, j)];
			v[IndexV(nx + 1, j)] = -v[IndexV(nx, j)];
		}

		//.....圧力方程式をガウスザイデル法で解く.....
		for (int k = 0; k < nk; k++) {
			double dpmax = 0.0;
			//境界の値を計算
			for (int i = 1; i < nx + 1; i++) {
				p[IndexP(i, 0)] = p[IndexP(i, 1)];
			}
			for (int i = 1; i < nx + 1; i++) {
				p[IndexP(i, ny + 1)] = p[IndexP(i, ny)];
			}
			for (int j = 1; j < ny + 1; j++) {
				p[IndexP(0, j)] = p[IndexP(1, j)];
			}
			for (int j = 1; j < ny + 1; j++) {
				p[IndexP(nx + 1, j)] = p[IndexP(nx, j)];
			}
			//内部点の計算
			for (int i = 1; i < nx + 1; i++) {
				for (int j = 1; j < ny + 1; j++) {
					if (i == 1 && j == 1) {
						p[IndexP(i, j)] = 0.0;
					} else {
						double dudx = (u[IndexU(i, j)] - u[IndexU(i - 1, j)])/dx, dvdy = (v[IndexV(i, j)] - v[IndexV(i, j - 1)])/dy, larged = dudx + dvdy, tmpp = p[IndexP(i, j)];
						p[IndexP(i, j)] = (1.0 - omega)*p[IndexP(i, j)] + omega*pow(dx*dy, 2.0)/(2.0*(pow(dx, 2.0) + pow(dy, 2.0)))*((p[IndexP(i + 1, j)] + p[IndexP(i - 1, j)])/pow(dx, 2.0) + (p[IndexP(i, j + 1)] + p[IndexP(i, j - 1)])/pow(dy, 2.0) - larged/dt + pow(larged, 2.0));
						if (dpmax < abs(tmpp - p[IndexP(i, j)])) {
							dpmax = abs(tmpp - p[IndexP(i, j)]);
						}
					}
				}
			}
			if (dpmax < eps) {
				cout << "\tConverged at k = " << k;
				break;
			}
		}

		//.....N-S方程式を解く.....
		//uについて
		for (int i = 1; i < nx; i++) {
			for (int j = 1; j < ny + 1; j++) {
				double ududx = (pow((u[IndexU(i + 1, j)] + u[IndexU(i, j)])/2.0, 2.0) - pow((u[IndexU(i, j)] + u[IndexU(i - 1, j)])/2.0, 2.0))/dx;
				double vdudy = ((u[IndexU(i, j + 1)] + u[IndexU(i, j)])*(v[IndexV(i + 1, j)] + v[IndexV(i, j)])/4.0 - (u[IndexU(i, j)] + u[IndexU(i, j - 1)])*(v[IndexV(i + 1, j - 1)] + v[IndexV(i, j - 1)])/4.0)/dy;

				double dpdx = (p[IndexP(i + 1, j)] - p[IndexP(i, j)])/dx;

				double d2udx2 = (u[IndexU(i + 1, j)] - 2.0*u[IndexU(i, j)] + u[IndexU(i - 1, j)])/pow(dx, 2.0);
				double d2udy2 = (u[IndexU(i, j + 1)] - 2.0*u[IndexU(i, j)] + u[IndexU(i, j - 1)])/pow(dy, 2.0);

				up1[IndexU(i, j)] = u[IndexU(i, j)] + (-ududx - vdudy - dpdx + (d2udx2 + d2udy2)/re)*dt;
			}
		}
		//vについて
		for (int i = 1; i < nx + 1; i++) {
			for (int j = 1; j < ny; j++) {
				double udvdx = ((u[IndexU(i, j + 1)] + u[IndexU(i, j)])*(v[IndexV(i + 1, j)] + v[IndexV(i, j)])/4.0 - (u[IndexU(i - 1, j + 1)] + u[IndexU(i - 1, j)])*(v[IndexV(i, j)] + v[IndexV(i - 1, j)])/4.0)/dx;
				double vdvdy = (pow((v[IndexV(i, j + 1)] + v[IndexV(i, j)])/2.0, 2.0) - pow((v[IndexV(i, j)] + v[IndexV(i, j - 1)])/2.0, 2.0))/dy;
					
				double dpdy = (p[IndexP(i, j + 1)] - p[IndexP(i, j)])/dy;

				double d2vdx2 = (v[IndexV(i + 1, j)] - 2.0*v[IndexV(i, j)] + v[IndexV(i - 1, j)])/pow(dx, 2.0);
				double d2vdy2 = (v[IndexV(i, j + 1)] - 2.0*v[IndexV(i, j)] + v[IndexV(i, j - 1)])/pow(dy, 2.0);

				vp1[IndexV(i, j)] = v[IndexV(i, j)] + (-udvdx - vdvdy - dpdy + (d2vdx2 + d2vdy2)/re)*dt;
			}
		}

		//速度の更新
		for (int i = 1; i < nx; i++) {
			for (int j = 1; j < ny + 1; j++) {
				u[IndexU(i, j)] = up1[IndexU(i, j)];
			}
		}
		for (int i = 1; i < nx + 1; i++) {
			for (int j = 1; j < ny; j++) {
				v[IndexV(i, j)] = vp1[IndexV(i, j)];
			}
		}
	}

	//----------結果の出力----------
	std::ofstream fout("result/mac.vtk");
    
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
    delete[] p, u, up1, v, vp1;

    return 0;
}