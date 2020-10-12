using System;
using System.Collections.Generic;
using System.IO;

namespace JPEG_LS
{
    public class Compress : Other
    {
        public void Compressing(byte[] data, int index, BinaryWriter stream, int height, int width)
        {
            // 1. Inicijalizacija promenljivih

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

            double[] arr_Errval = new double[512];
            double[] arr_MErrval = new double[512];

            for (int i = 0; i < height; i++)
            {
                for (int j = 0; j < width; j++)
                {
                    // Kontekstno modeliranje

                    double Ra, Rb, Rc, Rd, Ix = data[index];

                    Ra = (j == 0) ? 0 : data[index - 3];
                    Rb = (i == 0 || j == 0) ? 0 : data[index - (3 * width) - 3];
                    Rc = (i == 0) ? 0 : data[index - (3 * width)];
                    Rd = (i == 0 || j >= width - 1) ? 0 : data[index - (3 * width) + 3];

                    index += 3;


                    // A.1 Izračunavanje lokalnih gradijanata

                    double D1 = Rd - Rb,
                            D2 = Rb - Rc,
                            D3 = Rc - Ra;

                    // A.2 Selektovanje moda rada Run/Regular

                    if (Math.Abs(D1) <= NEAR && Math.Abs(D2) <= NEAR && Math.Abs(D3) <= NEAR)
                    {
                        #region goto RunModeProcessing

                        // RUN mod

                        //A.15 - Određivanje dužine trčanja

                        double RUNval = Ra,
                                RUNcnt = 0,
                                Rx = 0,
                                EOLine = 0;

                        while (Math.Abs(Ix - RUNval) <= NEAR)
                        {
                            ++RUNcnt;
                            Rx = RUNval;

                            if (j >= width - 2)
                            {
                                EOLine = 1;
                                break;
                            }
                            else
                            {
                                // GetNextSample () 
                                j++;
                                Ra = (j == 0) ? 0 : data[index - 3];
                                Rb = (i == 0 || j == 0) ? 0 : data[index - (3 * width) - 3];
                                Rc = (i == 0) ? 0 : data[index - (3 * width)];
                                Rd = (i == 0 || j >= width - 1) ? 0 : data[index - (3 * width) + 3];
                                Ix = data[index];
                                index += 3;
                            }
                        }

                        // A.16  Kodiranje dužine trčanja

                        int rm;

                        //A.16.1 Kodiranje dužine koja je rg 

                        while (RUNcnt >= (rm = J[RUNIndex]))
                        {
                            Write(stream, true);
                            RUNcnt = RUNcnt - rm;
                            if (RUNIndex < 31)
                            {
                                ++RUNIndex;
                            }
                        }

                        // A.16.2 Kodiranje dužine manje od rg 

                        if (Math.Abs(Ix - RUNval) > NEAR)
                        {
                            Write(stream, false);

                            int bits = J[RUNIndex],
                                    value = (int)RUNcnt;

                            while (Convert.ToBoolean(bits--))
                            {
                                bool bit = Convert.ToBoolean((value >> bits) & 1);

                                Write(stream, bit);
                            }

                            if (RUNIndex > 0)
                            {
                                --RUNIndex;
                            }
                        }
                        else if (RUNcnt > 0)
                        {
                            Write(stream, true);
                        }

                        // A.17  Prekid koji nije nastao krajem linije

                        if (EOLine != 1)
                        {
                            // A.17.1 Izračunavanje indeksa RItype koji definiše kontekst:

                            double RItype;

                            if (Math.Abs(Ra - Rb) <= NEAR)
                            {
                                RItype = 1;
                            }
                            else
                            {
                                RItype = 0;
                            }

                            //A.17.2   Izračunavanje greške predviđanja

                            double Px;

                            if (RItype == 1)
                            {
                                Px = Ra;
                            }
                            else
                            {
                                Px = Rb;
                            }

                            double Errval = Ix - Px, SIGN;

                            // A.17.3 Errval će biti kvantizovan i Rx izračunato, zatim će se greška smanjiti upotrebom promenljive RANGE

                            if ((RItype == 0) && (Ra > Rb))
                            {
                                Errval = -Errval;
                                SIGN = -1;
                            }
                            else
                            {
                                SIGN = 1;
                            }

                            if (NEAR > 0)
                            {
                                if (Errval > 0)
                                {
                                    Errval = (Errval + NEAR) / (2 * NEAR + 1);
                                }
                                else
                                {
                                    Errval = -((NEAR - Errval) / (2 * NEAR + 1));
                                }

                                Rx = Px + SIGN * Errval * (2 * NEAR + 1);
                            }
                            else
                            {
                                Rx = Ix;
                            }

                            if (Errval < 0)
                            {
                                Errval = Errval + RANGE;
                            }

                            if (Errval >= (RANGE + 1) / 2)
                            {
                                Errval = Errval - RANGE;
                            }

                            // A.17.4  Izračunavanje pomoćne promenljive TEMP

                            int Q = (int)(RItype + 365);

                            double TEMP;

                            if (RItype == 0)
                            {
                                TEMP = A[365];
                            }
                            else
                            {
                                TEMP = A[366] + ((int)N[366] >> 1);
                            }

                            double k = DetermineGolombParameter(N[Q], TEMP), EMErrval, map;

                            // A.17.5  Izračunavanje mape za mapiranje Errval:

                            if ((k == 0) && (Errval > 0) && (2 * Nn[Q - 365] < N[Q]))
                            {
                                map = 1;
                            }
                            else if ((Errval < 0) && (2 * Nn[Q - 365] >= N[Q]))
                            {
                                map = 1;
                            }
                            else if ((Errval < 0) && (k != 0))
                            {
                                map = 1;
                            }
                            else
                            {
                                map = 0;
                            }

                            EMErrval = 2 * Math.Abs(Errval) - RItype - map;

                            EncodeGolomb(k, LIMIT - J[RUNIndex] - 1, qbpp, EMErrval, ref stream);

                            // A.18 Ažuriranje promanljivih prekida: 

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

                        #endregion
                    }
                    else
                    {
                        #region goto RegularModeProcessing
                        // REGULAR mod

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

                        int Q = FunctionMappingVector(Q1, Q2, Q3);




                        // A.5 Izračunavanje vrednosti prediđanja Px: 
                        double Px = EdgeDetectingPredictor( Ra,  Rb, Rc);

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

                        // A.7 Izračunavanje greške predviđanja

                        double MErrval,
                                Errval = Ix - Px;

                        arr_Errval[(int)Errval + 255]++;

                        if (SIGN == -1)
                        {
                            Errval = -Errval;
                        }

                        // A.8  Kvantizacija greške: 

                        if (Errval > 0)
                        {
                            Errval = (Errval + NEAR) / (2 * NEAR + 1);
                        }
                        else
                        {
                            Errval = -(NEAR - Errval) / (2 * NEAR + 1);
                        }

                        double Rx = Px + SIGN * Errval * (2 * NEAR + 1);

                        Clipping(ref Rx, MAXVAL);

                        // A.9 Smanjenje greške: 

                        if (Errval < 0)
                        {
                            Errval = Errval + RANGE;
                        }

                        if (Errval >= ((RANGE + 1) / 2))
                        {
                            Errval = Errval - RANGE;
                        }

                        // A.10  Izračunavanje Golombove promenljive k:

                        double k = DetermineGolombParameter(N[Q], A[Q]);

                        // A.11  Mapiranje greške :

                        if ((NEAR == 0) && (k == 0) && (2 * B[Q] <= -N[Q]))
                        {
                            if (Errval >= 0)
                            {
                                MErrval = 2 * Errval + 1;
                            }
                            else
                            {
                                MErrval = -2 * (Errval + 1);
                            }
                        }
                        else
                        {
                            if (Errval >= 0)
                            {
                                MErrval = 2 * Errval;
                            }
                            else
                            {
                                MErrval = -2 * Errval - 1;
                            }
                        }

                        arr_MErrval[(int)MErrval + 255]++;

                        // A.12  Kodiranje mapirane vrednosti MErrval Golombovim kodiranjem: 
                        EncodeGolomb(k, LIMIT, qbpp, MErrval, ref stream);

                        //A.13 i A.14 Ažuriranje promenljivih:
                        //A.13

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
        }

        public void Write(BinaryWriter stream)
        {
            stream.Write(buffer);
            counter = 8;
        }
    }
}