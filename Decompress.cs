using System;
using System.Collections.Generic;

namespace JPEG_LS
{
    public class Decompress : Other
    {
        public double[,] Decompressing(List<bool> data, ref int index, int height, int width)
        {
            // Inicijalizacija

            double RANGE = Math.Ceiling((MAXVAL + 2 * NEAR) / (2 * NEAR + 1)) + 1;
            double qbpp = Math.Floor(Math.Log(RANGE, 2));
            double bpp = Math.Max(2, Math.Floor(Math.Log(MAXVAL + 1, 2)));
            double LIMIT = 2 * (bpp + Math.Max(8, bpp));

            
            double[] N = new double[nContexts + 2];
            double[] A = new double[nContexts + 2];
            double[] B = new double[nContexts];
            double[] C = new double[nContexts];

            for (int i = 0; i < nContexts + 2; i++)
            {
                N[i] = 1;
                A[i] = Math.Max(2, Math.Ceiling((RANGE + Math.Pow(2, 5)) / Math.Pow(2, 6)));
            }

            
            int[] J = { 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 9, 10, 11, 12, 13, 14, 15 };
            int RUNIndex = 0;

            
            double[] Nn = new double[2];

            double[,] component = new double[height, width];

            for (int i = 0; i < height; i++)
            {
                for (int j = 0; j < width; j++)
                {
                    double Ra, Rb, Rc, Rd;

                    Ra = (j == 0) ? 0 : component[i, j - 1];
                    Rb = (i == 0 || j == 0) ? 0 : component[i - 1, j - 1];
                    Rc = (i == 0) ? 0 : component[i - 1, j];
                    Rd = (i == 0 || j >= width - 1) ? 0 : component[i - 1, j + 1];

                    double D1 = Rd - Rb,
                            D2 = Rb - Rc,
                            D3 = Rc - Ra;

                    if (Math.Abs(D1) <= NEAR && Math.Abs(D2) <= NEAR && Math.Abs(D3) <= NEAR)
                    {
                        while (true)
                        {
                            int EOLine = 0;

                            if (data[index])
                            {
                                int rm = J[RUNIndex];

                                while (Convert.ToBoolean(rm--) && j < width - 1)
                                {
                                    component[i, j] = Ra;

                                    ++j;
                                }


                                if (rm == -1 && RUNIndex < 31)
                                {
                                    ++RUNIndex;
                                }

                                if (j >= width - 1)
                                {
                                    ++index;
                                    --j;
                                    break;
                                }

                            }
                            else
                            {
                                int bits = J[RUNIndex], value = 0, bit;

                                if (!data[index])
                                {
                                    index++;

                                    while (Convert.ToBoolean(bits--))
                                    {
                                        bit = Convert.ToInt32(data[index]);
                                        index++;

                                        value += bit << bits;
                                    }

                                    while (Convert.ToBoolean(value--))
                                    {
                                        if (j == width - 1)
                                        {
                                            break;
                                        }
                                        component[i, j] = Ra;
                                        ++j;
                                    }

                                    if (RUNIndex > 0)
                                    {
                                        --RUNIndex;
                                    }
                                }

                                if (EOLine != 1)
                                {
                                    // A.17.1 Odredjivanje parametra RIype
                                    double RItype;

                                    Ra = (j == 0) ? 0 : component[i, j - 1];
                                    Rb = (i == 0 || j == 0) ? 0 : component[i - 1, j - 1];

                                    if (Math.Abs(Ra - Rb) <= NEAR)
                                    {
                                        RItype = 1;
                                    }
                                    else
                                    {
                                        RItype = 0;
                                    }

                                    // A.17.2 Izračunavanje greške predviđanja: 
                                    double Px;

                                    if (RItype == 1)
                                    {
                                        Px = Ra;
                                    }
                                    else
                                    {
                                        Px = Rb;
                                    }

                                    double SIGN;

                                    if ((RItype == 0) && (Ra > Rb))
                                    {
                                        SIGN = -1;
                                    }
                                    else
                                    {
                                        SIGN = 1;
                                    }

                                    double TEMP;

                                    if (RItype == 0)
                                    {
                                        TEMP = A[365];
                                    }
                                    else
                                    {
                                        TEMP = A[366] + ((int)N[366] >> 1);
                                    }

                                    int Q = (int)(RItype + 365);

                                    double k = DetermineGolombParameter(N[Q], TEMP), Errval, EMErrval = 0;

                                    DecodeGolomb(k, LIMIT - J[RUNIndex] - 1, qbpp, ref EMErrval, data, ref index);


                                    double tEMErrval = EMErrval + RItype;

                                    if (tEMErrval == 0)
                                    {
                                        Errval = 0;
                                    }
                                    else if (k == 0)
                                    {
                                        if (2 * Nn[Q - 365] < N[Q])
                                        {
                                            if (tEMErrval % 2 == 0)
                                            {
                                                Errval = -((int)tEMErrval >> 1);        
                                            }
                                            else
                                            {
                                                Errval = ((int)tEMErrval + 1) >> 1;      
                                            }
                                        }
                                        else
                                        {   
                                            if (tEMErrval % 2 == 0)
                                            {
                                                Errval = (int)tEMErrval >> 1;            
                                            }
                                            else
                                            {
                                                Errval = -(((int)tEMErrval + 1) >> 1);  
                                            }
                                        }
                                    }
                                    else
                                    {
                                        if (tEMErrval % 2 == 0)
                                        {
                                            Errval = (int)tEMErrval >> 1;                
                                        }
                                        else
                                        {
                                            Errval = -(((int)tEMErrval + 1) >> 1);      
                                        }
                                    }

                                    if (SIGN < 0) Errval = -Errval;

                                    double Rx = Px + Errval;

                                    if (Rx < -NEAR)
                                    {
                                        Rx += RANGE * (2 * NEAR + 1);
                                    }
                                    else if (Rx > MAXVAL + NEAR)
                                    {
                                        Rx -= RANGE * (2 * NEAR + 1);
                                    }

                                    Clipping(ref Rx);

                                    component[i, j] = Rx;

                                    if (Errval < 0)
                                    {
                                        ++Nn[Q - 365];
                                    }

                                    A[Q] += ((int)(EMErrval + 1 - RItype) >> 1);

                                    if (N[Q] == RESET)
                                    {
                                        A[Q] = (int)A[Q] >> 1;
                                        N[Q] = (int)N[Q] >> 1;
                                        Nn[Q - 365] = (int)Nn[Q - 365] >> 1;
                                    }

                                    ++N[Q];
                                }

                                break;
                            }
                            index++;
                        }
                    }
                    else
                    {
                        #region goto RegularModeProcessing

                        // A.3 Kvantizacija lokalnih gradijenata (D1,D2,D3): 
                        int Q1 = LocalGradientQuantization(D1),
                            Q2 = LocalGradientQuantization(D2),
                            Q3 = LocalGradientQuantization(D3);

                        int SIGN;

                        // A.4 Spajanje kvantizovanih gradijenata:
                        if (Q1 < 0 || (Q1 == 0 && Q2 < 0) || (Q1 == 0 && Q2 == 0 && Q3 < 0))
                        {
                            Q1 = -Q1;
                            Q2 = -Q2;
                            Q3 = -Q3;
                            SIGN = -1;
                        }
                        else
                        {
                            SIGN = 1;
                        }

                        int Q;

                        if (Q1 == 0)
                        {
                            if (Q2 == 0)
                            {
                                Q = 360 + Q3;       
                            }
                            else
                            {   
                                Q = 324 + (Q2 - 1) * 9 + (Q3 + 4);  
                            }
                        }
                        else
                        {       
                            Q = (Q1 - 1) * 81 + (Q2 + 4) * 9 + (Q3 + 4);    
                        }

                        double Px;

                        // A.5 Izračunavanje vrednosti prediđanja Px: 

                        if (Rc >= Math.Max(Ra, Rb))
                        {
                            Px = Math.Min(Ra, Rb);
                        }
                        else if (Rc <= Math.Min(Ra, Rb))
                        {
                            Px = Math.Max(Ra, Rb);
                        }
                        else
                        {
                            Px = Ra + Rb - Rc;
                        }

                        // A.6 Korekcija predviđanja Px:

                        if (SIGN == 1)
                        {
                            Px = Px + C[Q];
                        }
                        else
                        {
                            Px = Px - C[Q];
                        }

                        Clipping(ref Px, MAXVAL);

                        double k = DetermineGolombParameter(N[Q], A[Q]);

                        double MErrval = 0,
                                Errval;

                        DecodeGolomb(k, LIMIT, qbpp, ref MErrval, data, ref index);

                        // A.11  Mapiranje greške :

                        if ((NEAR == 0) && (k == 0) && (2 * B[Q] <= -N[Q]))
                        {
                            if (MErrval % 2 != 0)
                            {
                                Errval = (MErrval - 1) / 2;
                            }
                            else
                            {
                                Errval = -(MErrval / 2) - 1;
                            }
                        }
                        else
                        {
                            if ((MErrval % 2) == 0)
                            {
                                Errval = MErrval / 2;
                            }
                            else
                            {
                                Errval = -(MErrval + 1) / 2;
                            }
                        }

                        double Rx = Px + SIGN * Errval;

                        if (Rx < -NEAR)
                        {
                            Rx += RANGE * (2 * NEAR + 1);
                        }
                        else if (Rx > MAXVAL + NEAR)
                        {
                            Rx -= RANGE * (2 * NEAR + 1);
                        }

                        Clipping(ref Rx);

                        component[i, j] = Rx;

                        // Ayurirranje promenljivih

                        B[Q] += Errval * (2 * NEAR + 1);
                        A[Q] += Math.Abs(Errval);



                        if (N[Q] == RESET)
                        {
                            A[Q] = (int)A[Q] >> 1;

                            if (B[Q] >= 0)
                            {
                                B[Q] = ((int)B[Q] >> 1);
                            }
                            else
                            {
                                B[Q] = -((1 - (int)B[Q]) >> 1);
                            }

                            N[Q] = (int)N[Q] >> 1;
                        }

                        ++N[Q];

                        // A.14

                        if (B[Q] <= -N[Q])
                        {
                            B[Q] += N[Q];

                            if (C[Q] > MIN_C)
                            {
                                --C[Q];
                            }

                            if (B[Q] <= -N[Q])
                            {
                                B[Q] = -N[Q] + 1;
                            }
                        }
                        else if (B[Q] > 0)
                        {
                            B[Q] -= N[Q];

                            if (C[Q] < MAX_C)
                            {
                                ++C[Q];
                            }

                            if (B[Q] > 0)
                            {
                                B[Q] = 0;
                            }
                        }
                        #endregion
                    }
                }
            }

            return component;
        }

    }
}