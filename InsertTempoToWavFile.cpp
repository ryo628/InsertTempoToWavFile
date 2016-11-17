// kadai 0707
// 拍への純音挿入プログラム

#define _CRT_SECURE_NO_DEPRECATE

#include<stdlib.h>
#include<stdio.h>
#include<string.h>

#include"../function_wave_fft.cpp"

#define ARG 3					// 入力時の文字数(例：a.exe input.wav output.wav)
#define TEMP_NUM 3				// 退避用の配列の数
#define N 8192					// 周波数解析に使うサンプル数(128, 256, 512, 1024, 2048, 4096, 8192, 16384, 32768)
#define ANALYSIS_START_SEC 1	// 解析を開始する秒数
#define ANALYSIS_END_SEC 60		// 解析を終了する秒数

int main(int argc, char *argv[])
{
	//エラー処理
	if(argc != ARG){fprintf(stderr, "Usage: program <InputWaveFile>\n");exit(1);}

	/* 変数の宣言 */
	// カウント用
	int i, j, n, point;

	// 周波数解析用
	int ip[NMAXSQRT + 2];
	double a[(4*N) + 1], w[(4*N) * 5 / 4];

	// Hz to Nへの変換係数用
	double aHz;

	// 解析を開始する時間及び解析を終了する時間の計算用
	int s_sec,e_sec;

	// 退避用
	double temp[TEMP_NUM];

	// 入力波形を構造体へ格納
	// 構造体を宣言
	Sound *snd;

	//入力(エラー処理)
	if((snd = Read_Wave(argv[1])) == NULL){exit(1);}

	//構造体Sound型内の変数「snd->****」の説明
	//snd->channelnum：入力されたwaveファイルのチャンネル数			例）ステレオならsnd->channelnum=2，モノラルならsnd->channelnum=1
	//snd->bit_per_sample：入力されたwaveファイルの量子化bit数		例）16bitならbit_per_sample=16
	//snd->samplinrate： 入力されたwaveファイルのサンプリングレート		例）サンプリングレート44100Hzのwaveファイルを読み込んだ場合snd->samplingrate=44100
	//snd->datanum：入力されたwaveファイルのデータ数				例）サンプリングレート44100Hz，2秒のwaveファイルを読み込んだ場合snd->datanum=88200

	// 波形のデータを格納する配列を宣言
	double *raw_wave = (double*)calloc(snd->datanum, sizeof(double));
	double *modified_wave = (double*)calloc(snd->datanum*10, sizeof(double));
	double *out_wave = (double*)calloc(snd->datanum, sizeof(double));

	// モノラルだった場合の処理（構造体からそのまま配列へ振幅値を格納）
	if(snd->channelnum == 1 && snd->bit_per_sample == 16)
	{
		for(i=0 ; i<snd->datanum ; i++) raw_wave[i] = snd[0].monaural16[i];
	}
	// ステレオだった場合の処理（構造体の左右のチャンネルを合成して配列へ振幅値を格納）
	else if(snd->channelnum == 2 && snd->bit_per_sample == 16)
	{
		for(i=0 ; i<snd->datanum ; i++) raw_wave[i] = (snd[0].stereo16[i].l/2) + (snd[0].stereo16[i].r/2);
	}
	// ステレオでもモノラルでもない場合のエラー処理
	else {printf("error : This file is unknown file type.");exit(1);}

	// Hz to n への変換係数の計算（Hz/aHz = i，i*aHz = Hz）
	aHz = (double) snd->samplingrate/N;


	//解析を開始する時刻と終了する時刻の計算
	s_sec = snd->samplingrate * ANALYSIS_START_SEC;
	e_sec = snd->samplingrate * ANALYSIS_END_SEC;

	//パワースペクトルを格納する配列の宣言
	double *power_01 = (double*)calloc(N/2, sizeof(double)); //出力用

	//音楽音響信号のテンポ推定
	//STEP1 : 高周波数成分の除去（ローパスフィルタに通す）
	for(point=0 ; point<snd->datanum-N ; point=point+(N/2))
	{
		//波形の入力＆窓関数をかける
		for(i=0 ; i<N ; i++){
			a[2*i]   = triangular_window(i,N) * raw_wave[point + i];	//フーリエ変換への入力信号の実部
			a[2*i+1] = 0.0;												//フーリエ変換への入力信号の虚部
		}

		//Fourier 変換（入力信号：入力した音響データ）
		ip[0] = 0;
		cdft(N*2, 1, a, ip, w);

		//ローパスフィルタに，フーリエ変換の結果を通す
		//パワースペクトルの左側に通す場合
		for(i=0 ; i<N/2-1 ; i++){
			if(110/aHz < i){
				a[2*i] = 0.0;
				a[2*i+1] = 0.0;
			}
		}
		//パワースペクトルの右側に通す場合
		for(i=N-1 ; i>N/2 ; i--){
			if((N-1)-(110/aHz) > i){
				a[2*i] = 0.0;
				a[2*i+1] = 0.0;
			}
		}

		//逆フーリエ変換
		ip[0] = 0;
		cdft(N*2, -1, a, ip, w);
		for (j=0 ; j<=2*N-1 ; j++) {
			a[j] *= 1.0 / N;
		}

		//ローパスフィルタに通した波形を別の配列に格納
		for(i=0 ; i<N ; i++){
			if( point+i <= snd->datanum){
				modified_wave[point+i] += a[2*i];
			}
		}

	}

	//STEP2 : ダウンサンプリングとパワー関数の算出
	//配列の初期化
	for(point=0 ; point<snd->datanum ; point++){ out_wave[point] = raw_wave[point]; raw_wave[point] = 0.0;}

	//ダウンサンプリング
	j=0;
	for(point=0 ; point<snd->datanum-200 ; point=point+200){
		raw_wave[j] = modified_wave[point];
		j++;
	}

	//パワー関数を算出
	for(point=0 ; point<j ; point++){
		raw_wave[point] *= raw_wave[point];
	}

	//STEP3 : テンポ推定とスペクトルの出力
	//ダウンサンプリング後のパワー関数の入力＆窓関数をかける
	for(i=0 ; i<N ; i++){
		a[2*i]   = triangular_window(i,N) * raw_wave[i];	//フーリエ変換への入力信号の実部
		a[2*i+1] = 0.0;										//フーリエ変換への入力信号の虚部
	}

	//Fourier 変換（入力信号：パワー関数）
	ip[0] = 0;cdft(N*2, 1, a, ip, w);

	//パワースペクトルの算出
	for(i=0 ; i<N/2 ; i++){
		power_01[i] = (a[2*i]*a[2*i])+(a[2*i+1]*a[2*i+1]);
	}

	//配列のインデックス（添え字）から周波数への変換係数の再計算（Hz/aHz = i，i*aHz = Hz）※ダウンサンプリングしたため
	aHz = (double)(snd->samplingrate/200)/N;

	// テンポ推定(maxに周波数代入)
	double tempo = 0.0;		// bpm
	int tempo_n = 10;		// 添え字
	// パワースペクトル最大値探索
	for(i=10; i<N/50 ; i++){ if( power_01[i] > power_01[tempo_n] ) tempo_n = i; }
	// テンポ推測
	for(i=10 ; i<N/50 ; i++)
	{
		// 拍子の調整
		// 4beat 8beatみたいなの判定のために適当な閾値を最初に超えたtempoを取り出す
		if( power_01[i] > power_01[tempo_n]/2 )
		{
			tempo = (i*aHz)*60;
			break;
		}
	}
	printf("tempo : %f\n", tempo); // debug用

	// パワー関数計算
	double *pow_wave = (double*)calloc(snd->datanum, sizeof(double));
	for( int i = 0; i < snd->datanum; i++ ) pow_wave[i] = (out_wave[i] * out_wave[i]);

	// パルス列 & 相互相関関数
	// ( パルス列必要ないから消した )
	int tempo_l = (int)( snd->samplingrate * 60 / tempo ); 			// テンポの周期のサンプル個数計算
	double *func_T = (double*)calloc( tempo_l, sizeof(double) );	// 相互相関関数
	printf("tempo len : %d\n", tempo_l); // debug用

	// 拍子の位相推測 & T最大計算
	int max_T = 0;
	for( int m = 0; m < tempo_l; m++ )
	{
		func_T[m] = 0.0;

		// 位相ずらしてtempo毎のスペクトル和算
		for( int n = m; n < snd->datanum; n += tempo_l ) func_T[m] += pow_wave[n];

		func_T[m] = func_T[m] / snd->datanum;

		if( func_T[m] > func_T[max_T] ) max_T = m;
	}
	printf("delta : %d\n", max_T); // debug用

	// メモリ解放( 重そうだから )
	free(pow_wave);
	free(func_T);

	// 純音挿入
	int pt = 600;									// 純音の周波数
	int teper = 100;								// テーパーの長さ
	int pt_l = 20 * snd->samplingrate / pt;			// 純音の長さ(20周期)
	double rad = ( 2 * PI ) / snd->samplingrate;	// rad
	// 曲全体ループ
	for( int i = max_T; i < snd->datanum; i += tempo_l )
	{
		// 1テンポ音ループ
		for( int j = -pt_l/2; j < pt_l/2; j++ )
		{
			// エラー処理
			if( ( 0 < i+j ) && ( i+j < snd->datanum ) )
			{
				// テーパー
				if( pt_l/2+j < teper )		out_wave[i+j] += 20000 * sin( rad * (j+pt_l/2) * pt ) *(pt_l/2+j)/teper;
				else if( pt_l/2-j < teper )	out_wave[i+j] += 20000 * sin( rad * (j+pt_l/2) * pt ) *(pt_l/2-j)/teper;
				else						out_wave[i+j] += 20000 * sin( rad * (j+pt_l/2) * pt );
			}
		}
	}

	// 正規化
	temp[0] = 0.0;
	//最大値抽出
	for( int i = 0; i < snd->datanum; i++ )
	{
		if( temp[0] < out_wave[i] ) temp[0] = out_wave[i];
	}
	//正規化
	for( int i = 0; i < snd->datanum; i++ )
	{
		out_wave[i] = 20000 * ( out_wave[i] / temp[0] );
	}

	// out
	for(point=snd->datanum ; point>0 ; point--){if(out_wave[point] != 0){break;}}
	temp[0] = (double)point/snd->samplingrate;

	//出力用の波形情報の宣言
	unsigned short channelnum     = 1;
	unsigned long samplingrate    = snd->samplingrate;
	unsigned short bit_per_sample = snd->bit_per_sample;
	double wave_sec = temp[0];

	unsigned long datasize = wave_sec * snd->samplingrate * 1 * (snd->bit_per_sample/8);
	Sound *snd_2;
	if((snd_2 = Create_Sound(channelnum, samplingrate, bit_per_sample, datasize)) == NULL){exit(1);}

	//波形の出力
	if(ARG == 3)
	{
		for( i=0 ; i<snd_2->datanum; i++ ) snd_2[0].monaural16[i] = out_wave[i];

		// 書き出し
		if(Write_Wave(argv[2], snd_2) != 0){exit(1);}
	}

	// 終了処理
	Free_Sound(snd);
	Free_Sound(snd_2);
	return 0;
}
