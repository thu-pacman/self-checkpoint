/* 
 * -- High Performance Computing Linpack Benchmark (HPL)                
 *    HPL - 2.2 - February 24, 2016                          
 *    Antoine P. Petitet                                                
 *    University of Tennessee, Knoxville                                
 *    Innovative Computing Laboratory                                 
 *    (C) Copyright 2000-2008 All Rights Reserved                       
 *                                                                      
 * -- Copyright notice and Licensing terms:                             
 *                                                                      
 * Redistribution  and  use in  source and binary forms, with or without
 * modification, are  permitted provided  that the following  conditions
 * are met:                                                             
 *                                                                      
 * 1. Redistributions  of  source  code  must retain the above copyright
 * notice, this list of conditions and the following disclaimer.        
 *                                                                      
 * 2. Redistributions in binary form must reproduce  the above copyright
 * notice, this list of conditions,  and the following disclaimer in the
 * documentation and/or other materials provided with the distribution. 
 *                                                                      
 * 3. All  advertising  materials  mentioning  features  or  use of this
 * software must display the following acknowledgement:                 
 * This  product  includes  software  developed  at  the  University  of
 * Tennessee, Knoxville, Innovative Computing Laboratory.             
 *                                                                      
 * 4. The name of the  University,  the name of the  Laboratory,  or the
 * names  of  its  contributors  may  not  be used to endorse or promote
 * products  derived   from   this  software  without  specific  written
 * permission.                                                          
 *                                                                      
 * -- Disclaimer:                                                       
 *                                                                      
 * THIS  SOFTWARE  IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES,  INCLUDING,  BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE UNIVERSITY
 * OR  CONTRIBUTORS  BE  LIABLE FOR ANY  DIRECT,  INDIRECT,  INCIDENTAL,
 * SPECIAL,  EXEMPLARY,  OR  CONSEQUENTIAL DAMAGES  (INCLUDING,  BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA OR PROFITS; OR BUSINESS INTERRUPTION)  HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT,  STRICT LIABILITY,  OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE. 
 * ---------------------------------------------------------------------
 */ 

/*
 * --------------------------------------
 * The SKT-HPL project modifies this file
 * ---------------------------------------
 */
#define LOGL if(GRID->iam == 0) printf("LINE: %d\n", __LINE__);

/*
 * Include files
 */
#include "hpl.h"

/*
 * For snapshot
 */
#include <stdio.h>
#include <string.h>
#include <sys/ipc.h>
#include <sys/shm.h>

#ifdef STDC_HEADERS
  void HPL_pdgesv0
(
 HPL_T_grid *                     GRID,
 HPL_T_palg *                     ALGO,
 HPL_T_pmat *                     A
 )
#else
  void HPL_pdgesv0
( GRID, ALGO, A )
  HPL_T_grid *                     GRID;
  HPL_T_palg *                     ALGO;
  HPL_T_pmat *                     A;
#endif
{
  /* 
   * Purpose
   * =======
   *
   * HPL_pdgesv0 factors a N+1-by-N matrix using LU factorization with row
   * partial pivoting.  The main algorithm  is the "right looking" variant
   * without look-ahead. The lower triangular factor is left unpivoted and
   * the pivots are not returned. The right hand side is the N+1 column of
   * the coefficient matrix.
   *
   * Arguments
   * =========
   *
   * GRID    (local input)                 HPL_T_grid *
   *         On entry,  GRID  points  to the data structure containing the
   *         process grid information.
   *
   * ALGO    (global input)                HPL_T_palg *
   *         On entry,  ALGO  points to  the data structure containing the
   *         algorithmic parameters.
   *
   * A       (local input/output)          HPL_T_pmat *
   *         On entry, A points to the data structure containing the local
   *         array information.
   *
   * ---------------------------------------------------------------------
   */ 
  /*
   * .. Local Variables ..
   */
  HPL_T_panel                * * panel = NULL;
  HPL_T_UPD_FUN              HPL_pdupdate;
  int                        N, j, jb, n, nb, tag=MSGID_BEGIN_FACT,
                             test=HPL_KEEP_TESTING;
  int                        myrow, mycol, nprow, npcol;
#ifdef HPL_PROGRESS_REPORT
  double start_time, time, gflops;
#endif
  /* ..
   * .. Executable Statements ..
   */
  if( ( N = A->n ) <= 0 ) return;

#ifdef HPL_PROGRESS_REPORT
  start_time = HPL_timer_walltime();
#endif
  int sum_rank = GRID_sum_rank;
  int sum_size = GRID_sum_size;
  if(( N = A->n ) <= 0 ) return;

  HPL_grid_info(GRID, &nprow, &npcol, &myrow, &mycol);
  size_t Acnt = (size_t)(A->ld)*(size_t)(A->nq);
  double* Bptr = &A->A[Acnt - A->ld];
  size_t Scnt = (Acnt-1)/(sum_size-1)+1; // strip size
  Acnt = Scnt * sum_size; // Acnt = ceil[Acnt/(T-1)] * T
  if (GRID->iam == 0) 
  {
    printf("Acnt = %d, Scnt = %d, A->ld = %d, A->nq = %d\n", Acnt, Scnt, A->ld, A->nq);
  }

  /*
   * Allocate space for snapshot
   */
  double bt0 = MPI_Wtime();
  char* shm = NULL;
  key_t key = GRID->iam + 0x6789; // shm key
  int shmid;
  // the shm should be build before xhpl launch
  size_t shmsz = (Acnt + 64)*sizeof(double);
  if ((shmid = shmget(key, shmsz, IPC_CREAT | 0666)) < 0) { perror("shmget"); exit(1); }
  if ((shm = shmat(shmid, NULL, 0)) == (char *) -1) { perror("shmat"); exit(1); }
  char* pValid = shm;
  double *pA = (double*)HPL_PTR(&shm[1], ALGO->align);
  double *Asum = &(A->A[Acnt-Scnt]);

  int blksz = (1<<20); //xM double = 128x MB of the whole node (16 procs)
  //double *Asum = (double*)malloc((Scnt+ALGO->align) * sizeof(double));
  //Asum = (double*)HPL_PTR(Asum, ALGO->align);

  int do_recover = 0;
  char* recover_envstr = getenv("RECOVER");
  if (recover_envstr != NULL && atoi(recover_envstr) > 0) do_recover = 1;
  int only_one_snapshot = 0;
  char* only_envstr = getenv("ONLY");
  if (only_envstr != NULL && atoi(only_envstr) > 0) only_one_snapshot = 1;
  int do_snapshot = 0;
  int snapshot_interval = -1; // j nbs between snapshot
  char* snapshot_envstr = getenv("SNAPSHOT");
  if (snapshot_envstr != NULL && (snapshot_interval = atoi(snapshot_envstr)) > 0) {
    do_snapshot = 1;
  }

  int RUN_TIME_LIMIT = atoi(getenv("RUN_TIME_LIMIT"));

  int has_lost = 0, lost_rank = -1;
  int Lpass = 0;
  double last_t, this_t;
  last_t = MPI_Wtime();

  // set b to zero
  if (mycol != 0) {
    memset(Bptr, 0, A->ld*sizeof(double));
  }
  double bt1 = MPI_Wtime();
  if (GRID->iam == 0)
  {
    printf("PREPARING WORK COSTS %8.2lf seconds\n", bt1-bt0);
    printf("ITVL = %d\n", snapshot_interval);
  }

  HPL_pdupdate = ALGO->upfun; nb = A->nb;
  /*
   * Allocate a panel list of length 1 - Allocate panel[0] resources
   */
  panel = (HPL_T_panel **)malloc( sizeof( HPL_T_panel * ) );
  if( panel == NULL )
  { HPL_pabort( __LINE__, "HPL_pdgesv0", "Memory allocation failed" ); }

  HPL_pdpanel_new( GRID, ALGO, N, N+1, Mmin( N, nb ), A, 0, 0, tag,
      &panel[0] );
  /*
   * Loop over the columns of A
   */
  for( j = 0; j < N; j += nb )
  {
    //printf("LOOP No.%d, j=%d\n", j/nb, j);
    /*
     * Recover from the snapshot
     * A, tag, j,
     * start_time, dprint
     */ 
    if (do_recover) {
      double rt0 = MPI_Wtime();
      do_recover = 0;

      if (*pValid == '*') // good checkpoint
      {
        memcpy(A->A, pA, (Acnt-Scnt)*sizeof(double));
        memcpy(Asum, &pA[Acnt-Scnt], Scnt*sizeof(double));
        j = *(int*)&pA[Acnt];
        if (GRID->iam == 0) printf("Get j from snapshot, j = %d\n", j);
        //printf("!!!!!!   Performing Recover from memory id:%d, j = %d\n", shmid, j);
      }
#if 1
      else if (*pValid == '#') // good workspace
      {
        j = *(int*)&pA[Acnt];
        if (GRID->iam == 0) printf("Get j from workspace, j = %d\n", j);
      }
#endif
      else {
        //printf("RANK:%d The shm doesnt contain a valid snapshot\n", GRID->iam);
        has_lost = 1;
        lost_rank = sum_rank;
      }
      double rt1 = MPI_Wtime();
      MPI_Allreduce(MPI_IN_PLACE, &has_lost, 1, MPI_INT, MPI_SUM, GRID_sum_comm);
      int max_lost_cnt = has_lost;
      MPI_Allreduce(MPI_IN_PLACE, &max_lost_cnt, 1, MPI_INT, MPI_MAX, GRID->all_comm);
      //if (has_lost > 0) printf("There is %d lost in my family\n", has_lost);
      if (max_lost_cnt > 1) {
        if (GRID->iam == 0) printf("THIS IS A FRESH RUN (OR TOO MANY LOST), START FROM THE BEGINNING\n");
        j = -nb;
      }
      else if (has_lost == 1) {
        MPI_Allreduce(MPI_IN_PLACE, &lost_rank, 1, MPI_INT, MPI_MAX, GRID_sum_comm);
        if (lost_rank == sum_rank) { // I'm the lost rank
          //printf("Try to recover single rank failure, sum_rank = %d (Global rank = %d)\n", lost_rank, GRID->iam); 
          memset(A->A, 0, (Acnt-Scnt) * sizeof(double));
          memset(Asum, 0, Scnt * sizeof(double));
          int strip;
          for (strip = 0; strip < sum_size; ++ strip) {
            if (strip == sum_rank) {
              int bk;
              for (bk = Scnt; bk > 0; bk -= blksz) {
                int Asz = (bk >= blksz)? blksz: bk;
                MPI_Reduce(MPI_IN_PLACE, &Asum[Scnt - bk], Asz, MPI_DOUBLE, MPI_SUM, lost_rank, GRID_sum_comm);
              }
            }
            else {
              int d_strip = (strip > sum_rank) ? strip-1 : strip;
              int bk;
              for (bk = Scnt; bk > 0; bk -= blksz) {
                int Asz = (bk >= blksz)? blksz: bk;
                MPI_Reduce(MPI_IN_PLACE, &A->A[d_strip*Scnt + Scnt - bk], Asz, MPI_DOUBLE, MPI_SUM, lost_rank, GRID_sum_comm);
              }
              int aidx;
            }
          }
          int aidx;
          for (aidx = 0; aidx < Acnt - Scnt; ++ aidx) {
            A->A[aidx] *= -1;
          }
        }
        else { // healthy ranks
          int aidx;
          for (aidx = 0; aidx < Scnt; ++ aidx) {
            Asum[aidx] *= -1;
          }
          int strip;
          for (strip = 0; strip < sum_size; ++ strip) {
            if (strip == sum_rank) {
              int bk;
              for (bk = Scnt; bk > 0; bk -= blksz) {
                int Asz = (bk >= blksz)? blksz: bk;
                MPI_Reduce(&Asum[Scnt - bk], NULL, Asz, MPI_DOUBLE, MPI_SUM, lost_rank, GRID_sum_comm);
              }
            }
            else {
              int d_strip = (strip > sum_rank) ? strip-1 : strip;
              int bk;
              for (bk = Scnt; bk > 0; bk -= blksz) {
                int Asz = (bk >= blksz)? blksz: bk;
                MPI_Reduce(&A->A[d_strip*Scnt + Scnt - bk], NULL, Asz, MPI_DOUBLE, MPI_SUM, lost_rank, GRID_sum_comm);
              }
              int aidx;
            }
          }
          for (aidx = 0; aidx < Scnt; ++ aidx) {
            Asum[aidx] *= -1;
          }
        }

        int scalar_sender = (lost_rank + 1) % sum_size;
        if (sum_rank == scalar_sender) {
          //printf("I am scalar_sender, send to %d\n", lost_rank);
          MPI_Send(&j, 1, MPI_INT, lost_rank, 0, GRID_sum_comm);
        }
        if (sum_rank == lost_rank) {
          MPI_Recv(&j, 1, MPI_INT, scalar_sender, 0, GRID_sum_comm, MPI_STATUS_IGNORE);
        }
        // rewrite the snapshot
        // if (sum_rank == lost_rank) 
        {
          //    printf("Single rank failure fixed\n");
          //printf("REWRITE THE LOST SNAPSHOT, J=%d\n", j);
          *pValid = '#';
          memcpy(pA, A->A, (Acnt-Scnt)*sizeof(double));
          memcpy(&pA[Acnt-Scnt], Asum, Scnt*sizeof(double));
          *(int*)&pA[Acnt] = j;
          *pValid = '*';
        }
      }
      // J in lost_rank will recover from other ranks
      j += nb; if (j >= N) break;
      double rt2 = MPI_Wtime();
      double rt3 = rt1-rt0;
      MPI_Allreduce(MPI_IN_PLACE, &rt1, 1, MPI_DOUBLE, MPI_MIN, GRID->all_comm);
      MPI_Allreduce(MPI_IN_PLACE, &rt2, 1, MPI_DOUBLE, MPI_MAX, GRID->all_comm);
      if (max_lost_cnt <= 1 && GRID->iam  == 0) 
        printf("THIS RECOVERY COSTS %8.2lf (load-data), %8.2lf (rebuild-snapshot), Total %.2lf sec\n", rt3, rt2-rt1, rt2-rt0);
    }

#ifdef ENDEARLY
    if ( dprint >= .04 )
    {
      A->info = j;
      return ;
    }
#endif
    n = N - j; jb = Mmin( n, nb );
#ifdef HPL_PROGRESS_REPORT
    /* if this is process 0,0 and not the first panel */
    if ( GRID->myrow == 0 && GRID->mycol == 0 && j > 0 ) 
    {
      time = HPL_timer_walltime() - start_time;
      gflops = 2.0*(N*(double)N*N - n*(double)n*n)/3.0/(time > 0.0 ? time : 1e-6)/1e9;
      HPL_fprintf( stdout, "Column=%09d Fraction=%4.1f%% Gflops=%9.3e\n", j, j*100.0/N, gflops);
    }
#endif
    /*
     * Release panel resources - re-initialize panel data structure
     */
    (void) HPL_pdpanel_free( panel[0] );
    HPL_pdpanel_init( GRID, ALGO, n, n+1, jb, A, j, j, tag, panel[0] );
    /*
     * Factor and broadcast current panel - update
     */
    HPL_pdfact(               panel[0] );
    (void) HPL_binit(         panel[0] );
    do
    { (void) HPL_bcast(       panel[0], &test ); }
    while( test != HPL_SUCCESS );
    (void) HPL_bwait(         panel[0] );
    HPL_pdupdate( NULL, NULL, panel[0], -1 );
    /*
     * Update message id for next factorization
     */
    tag = MNxtMgid( tag, MSGID_BEGIN_FACT, MSGID_END_FACT );
    /*
     * build a snapshot
     * A, tag, j,
     * start_time, dprint
     */
    this_t = MPI_Wtime();
    if (do_snapshot && j*100/N < 60 && j/nb%snapshot_interval == 0) {
      *pValid = '#';
      if (only_one_snapshot) do_snapshot = 0;
      // the leading Aoff cols didnt update since last checkpoint
      int Aoff = Lpass*A->ld;
      double *Anew = A->A, *Aold = pA;
      // use the increamental checkpoint

      double st0 = MPI_Wtime();
      // set the sum-zone to 0
      memset(Asum, 0, Scnt * sizeof(double));
      double st1 = MPI_Wtime();
      int strip;
      for (strip = Aoff/Scnt; strip < sum_size; ++ strip) {
        if (strip == sum_rank) {
          int bk;
          for (bk = Scnt; bk > 0; bk -= blksz) {
            int Asz = (bk >= blksz)? blksz: bk;
            MPI_Reduce(MPI_IN_PLACE, &Asum[Scnt - bk], Asz, MPI_DOUBLE, MPI_SUM, strip, GRID_sum_comm);
          }
        }
        else {
          int d_strip = (strip > sum_rank) ? strip-1 : strip;
          int bk;
          for (bk = Scnt; bk > 0; bk -= blksz) {
            int Asz = (bk >= blksz)? blksz: bk;
            MPI_Reduce(&Anew[d_strip*Scnt + Scnt - bk], NULL, Asz, MPI_DOUBLE, MPI_SUM, strip, GRID_sum_comm);
          }
        }
      }
      double st2 = MPI_Wtime();

      *(int*)&pA[Acnt] = j;
#ifdef ASYOUGO
      pA[Acnt+1] = start_time;
      pA[Acnt+2] = dprint;
#endif
      if (GRID->iam == 0) printf("partially overwrite chkpt data\n");
      sleep(2);
      memcpy(&Aold[Aoff], &Anew[Aoff], (Acnt-Scnt-Aoff)*sizeof(double));
      double st3 = MPI_Wtime();

      if (Aoff/Scnt <= sum_rank) 
        memcpy(&Aold[Acnt-Scnt], Asum, Scnt*sizeof(double));
      double st4 = MPI_Wtime();

      *pValid = '*';

      double st_ms = st1-st0;
      double st_rd = st2-st1;
      double st_cm = st3-st2;
      double st_cs = st4-st3;
      MPI_Allreduce(MPI_IN_PLACE, &st_ms, 1, MPI_DOUBLE, MPI_MAX, GRID->all_comm);
      MPI_Allreduce(MPI_IN_PLACE, &st_rd, 1, MPI_DOUBLE, MPI_MAX, GRID->all_comm);
      MPI_Allreduce(MPI_IN_PLACE, &st_cm, 1, MPI_DOUBLE, MPI_MAX, GRID->all_comm);
      MPI_Allreduce(MPI_IN_PLACE, &st_cs, 1, MPI_DOUBLE, MPI_MAX, GRID->all_comm);

      if (GRID->iam == 0) 
      {
        printf("SNAPSHOT %6d MB/rank, MEMSET %5.2lf sec, ENCODE %5.2lf sec, %4.2lf (cpy-mat), %4.2lf (cpy-sum), Aoff = %d\n", 
            ((Acnt-Aoff)*sizeof(double))>>20, st_ms, st_rd, st_cm, st_cs, Aoff);
      }

      Lpass = j/(nb*npcol) * nb;
      last_t = MPI_Wtime();
      //if (sum_rank == 0) printf("put j=%d into snapshot\n", *(int*)&pA[Acnt]);
    }

    if ( myrow==0 && mycol==0)
    {
      if (HPL_timer_walltime() - start_time > RUN_TIME_LIMIT)
      {
        double tt = HPL_timer_walltime() - start_time;
        double dtmp = (double) N;
        double dtmp1 = (double)(N-j);
        double mflops = 2.0*(dtmp*dtmp*dtmp-dtmp1*dtmp1*dtmp1)/3.0;
        mflops = mflops / (1000000.0*tt);

        printf("RUN_TIME(%d s) END, MFLOPS=%8.2f\n", RUN_TIME_LIMIT, mflops);
        fflush(NULL);

        MPI_Finalize();
        exit(0);
      }
    }
  }
  /*
   * Release panel resources and panel list
   */
  (void) HPL_pdpanel_disp( &panel[0] );

  if( panel ) free( panel );
  /*
   * End of HPL_pdgesv0
   */
}
