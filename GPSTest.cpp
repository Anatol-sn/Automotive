#include "stdafx.h"
//#include "resource.h"

#include <SetupAPI.h>
#include <devguid.h>
#pragma comment( lib, "Setupapi.lib")
#include "Tools.h"		
#include "GeoDeclination.h"
/*
MINMEA
 
1. Define MINMEA_INCLUDE_COMPAT in the build environment.
If your GPS receiver outputs very long sentences, consider increasing MINMEA_MAX_SENTENCE_LENGTH
*/

#include "minmea\minmea.h"
#define INDENT_SPACES "  "

LPCWSTR pszVID = L"VID_1546";  // VendorID for gps UBLOX

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
struct MyCoordinates {
	bool   m_fValid;
	double m_dLatitude;
	double m_dLongitude;
	double m_dHeight;

	MyCoordinates() : m_fValid(false)
	{}
};


//---------------------------------------------------------------------------
int MiNMEA( LPCSTR line, MyCoordinates& coor )
{
	switch (minmea_sentence_id(line, false)) {
	case MINMEA_SENTENCE_GBS: {
		struct minmea_sentence_gbs frame;
		if (minmea_parse_gbs(&frame, line)) {
			printf(INDENT_SPACES "*$xxGBS: errors %f|%f|%f in meters\n",
				minmea_tofloat( &frame.err_latitude ),
				minmea_tofloat( &frame.err_longitude ),
				minmea_tofloat( &frame.err_altitude ) );
		}
		else {
			printf(INDENT_SPACES "$xxRMC sentence is not parsed\n");
		}
	} break;

	case MINMEA_SENTENCE_RMC: {
		struct minmea_sentence_rmc frame;
		if (minmea_parse_rmc(&frame, line)) {
			printf(INDENT_SPACES "*$xxRMC: %s, coords %f|%f, date %d.%d.%d\n",
				frame.valid ? "ok" : "INVALID",
				minmea_tocoord(&frame.latitude),
				minmea_tocoord(&frame.longitude),
				frame.date.year, 
				frame.date.month,
				frame.date.day
				);
		}
		else {
			printf(INDENT_SPACES "$xxRMC sentence is not parsed\n");
		}
	} break;

	case MINMEA_SENTENCE_GGA: {
		struct minmea_sentence_gga frame;
		if (minmea_parse_gga(&frame, line)) {
			printf(INDENT_SPACES "*$xxGGA: coords %f|%f|%f%c, height %f%c, fix quality %d, sats %d, hdop %f\n", 
				minmea_tocoord(&frame.latitude),
				minmea_tocoord(&frame.longitude),
				minmea_tofloat(&frame.altitude), frame.altitude_units,
				minmea_tofloat(&frame.height), frame.height_units,
				frame.fix_quality,
				frame.satellites_tracked,
				minmea_tofloat( &frame.hdop )
//					minmea_tofloat( &frame.dgps_age ) 
			);
/*
qualities
0: Fix not valid
1: GPS fix
2: Differential GPS fix (DGNSS), SBAS, OmniSTAR VBS, Beacon, RTX in GVBS mode
3: Not applicable
4: RTK Fixed, xFill, RTX
5: RTK Float, OmniSTAR XP/HP, Location RTK
6: INS Dead reckoning 
*/
			coor.m_dHeight = minmea_tofloat(&frame.altitude);
			coor.m_dLatitude = minmea_tocoord(&frame.latitude);
			coor.m_dLongitude = minmea_tocoord(&frame.longitude);
			coor.m_fValid = true;
		}
		else {
			printf(INDENT_SPACES "$xxGGA sentence is not parsed\n");
		}
	} break;

	case MINMEA_SENTENCE_GST: {
		struct minmea_sentence_gst frame;
		if (minmea_parse_gst(&frame, line)) {
			printf(INDENT_SPACES "*$xxGST lat, long, alt error deviation %fm|%fm|%fm\n",
				minmea_tofloat(&frame.latitude_error_deviation),
				minmea_tofloat(&frame.longitude_error_deviation),
				minmea_tofloat(&frame.altitude_error_deviation));
		}
		else {
			printf(INDENT_SPACES "$xxGST sentence is not parsed\n");
		}
	} break;

	case MINMEA_SENTENCE_GSV: {
		struct minmea_sentence_gsv frame;
		if (minmea_parse_gsv(&frame, line)) {
			printf(INDENT_SPACES "$xxGSV: message %d of %d\n", frame.msg_nr, frame.total_msgs);
			printf(INDENT_SPACES "$xxGSV: satellites in view: %d\n", frame.total_sats);
			for (int i = 0; i < 4; i++)
				printf(INDENT_SPACES "$xxGSV: sat nr %d, elevation: %d, azimuth: %d, snr: %d dbm\n",
					frame.sats[i].nr,
					frame.sats[i].elevation,
					frame.sats[i].azimuth,
					frame.sats[i].snr);
		}
		else {
			printf(INDENT_SPACES "$xxGSV sentence is not parsed\n");
		}
	} break;

	case MINMEA_SENTENCE_VTG: {
		struct minmea_sentence_vtg frame;
		if (minmea_parse_vtg(&frame, line)) {
			printf(INDENT_SPACES "$xxVTG: true track degrees = %f\n",
				minmea_tofloat(&frame.true_track_degrees));
			printf(INDENT_SPACES "        magnetic track degrees = %f\n",
				minmea_tofloat(&frame.magnetic_track_degrees));
			printf(INDENT_SPACES "        speed knots = %f\n",
				minmea_tofloat(&frame.speed_knots));
			printf(INDENT_SPACES "        speed kph = %f\n",
				minmea_tofloat(&frame.speed_kph));
		}
		else {
			printf(INDENT_SPACES "$xxVTG sentence is not parsed\n");
		}
	} break;

	case MINMEA_SENTENCE_ZDA: {
		struct minmea_sentence_zda frame;
		if (minmea_parse_zda(&frame, line)) {
			printf(INDENT_SPACES "$xxZDA: %d:%d:%d %02d.%02d.%d UTC%+03d:%02d\n",
				frame.time.hours,
				frame.time.minutes,
				frame.time.seconds,
				frame.date.day,
				frame.date.month,
				frame.date.year,
				frame.hour_offset,
				frame.minute_offset);
		}
		else {
			printf(INDENT_SPACES "$xxZDA sentence is not parsed\n");
		}
	} break;

	case MINMEA_INVALID: {
		printf(INDENT_SPACES "$xxxxx sentence is not valid\n");
	} break;

	default: {
		printf(INDENT_SPACES "$xxxxx sentence is not parsed\n");
	} break;
	}

	return 0;
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

volatile bool g_fExit = false;

// note: it is recommended to catch WM_QUERYENDSESSION/WM_ENDSESSION, if windows (GUI) will be added to this process;
BOOL WINAPI ConsoleHandler( DWORD dwCtrlType )
{
	if ( CTRL_SHUTDOWN_EVENT == dwCtrlType )
	{
		wprintf_s(L"system is going to shut down, finishing..." );
		g_fExit = true;
		return TRUE;
	}
	else if ( CTRL_LOGOFF_EVENT == dwCtrlType )
	{
		wprintf_s(L"user is logging off, finishing..." );
		g_fExit = true;
		return TRUE;
	}
	else if ( CTRL_CLOSE_EVENT == dwCtrlType || CTRL_C_EVENT == dwCtrlType || CTRL_BREAK_EVENT == dwCtrlType )
	{
		wprintf_s(L"user asks to quit, finishing..." );
		g_fExit = true;
		return TRUE;
	}

	return FALSE;
}
//---------------------------------------------------------------------------

void cls(HANDLE hConsole)
{
	CONSOLE_SCREEN_BUFFER_INFO csbi;
	SMALL_RECT scrollRect;
	COORD scrollTarget;
	CHAR_INFO fill;

	// Get the number of character cells in the current buffer.
	if (!GetConsoleScreenBufferInfo(hConsole, &csbi))
	{
		return;
	}

	// Scroll the rectangle of the entire buffer.
	scrollRect.Left = 0;
	scrollRect.Top = 0;
	scrollRect.Right = csbi.dwSize.X;
	scrollRect.Bottom = csbi.dwSize.Y;

	// Scroll it upwards off the top of the buffer with a magnitude of the entire height.
	scrollTarget.X = 0;
	scrollTarget.Y = (SHORT)(0 - csbi.dwSize.Y);

	// Fill with empty spaces with the buffer's default text attribute.
	fill.Char.UnicodeChar = TEXT(' ');
	fill.Attributes = csbi.wAttributes;

	// Do the scroll
	ScrollConsoleScreenBuffer(hConsole, &scrollRect, NULL, scrollTarget, &fill);

	// Move the cursor to the top left corner too.
	csbi.dwCursorPosition.X = 0;
	csbi.dwCursorPosition.Y = 0;

	SetConsoleCursorPosition(hConsole, csbi.dwCursorPosition);
}
//---------------------------------------------------------------------------

void ConvertByteToWChar(const BYTE* lpByte, DWORD dwSize, wchar_t* wchrValue)
{
	DWORD i = 0;
	DWORD index = 0;
	for (i = 0; i < dwSize; i += 2)
	{
		wchrValue[index++] = lpByte[i];
	}
}
//---------------------------------------------------------------------------

HRESULT DoMain()
{
	CDeclinationProc Declination;

	WCHAR szRealPath[MAX_PATH] = { 0 };
	
	HDEVINFO hdi = NULL;
	hdi = SetupDiGetClassDevs((const GUID*)&GUID_DEVCLASS_PORTS, NULL, NULL, DIGCF_PRESENT);

	if (hdi != INVALID_HANDLE_VALUE)
	{
		SP_DEVINFO_DATA dd = { 0 };
		dd.cbSize = sizeof(dd);

		for (UINT idx = 0; SetupDiEnumDeviceInfo(hdi, idx, &dd); idx++)
		{
			WCHAR szInst[MAX_PATH] = { 0 };
			DWORD dws = 0;
			if (!SetupDiGetDeviceInstanceId(hdi, &dd, szInst, MAX_PATH, &dws))
				continue;

			if (!wcsstr(szInst, pszVID))
			{
				//msyslog(LOG_INFO, "VAR: spec instance %s not found in %s, check next", szSpecInstance, szInst);
				continue;
			}

			//msyslog(LOG_INFO, "VAR: spec instance %s found in %s", szSpecInstance, szInst);

			HKEY key = SetupDiOpenDevRegKey(hdi, &dd, DICS_FLAG_GLOBAL, 0, DIREG_DEV, KEY_QUERY_VALUE);
			if (key == INVALID_HANDLE_VALUE)
			{
				//msyslog(LOG_ERR, "VAR: failed to open aliased device regkey, error %u", GetLastError());
				continue;
			}

			BYTE bVal[2*MAX_PATH] = { 0 };
			WCHAR szVal[MAX_PATH] = { 0 };
			ULONG ulLen = MAX_PATH;
			DWORD dwType;

			LRESULT lr = RegQueryValueEx(key, L"PortName", NULL, &dwType, bVal, &ulLen);
			
			ConvertByteToWChar(bVal, ulLen, szVal);
			
			RegCloseKey(key);

			if (lr == ERROR_SUCCESS && (dwType == REG_SZ || dwType == REG_EXPAND_SZ))
			{
				wcscpy_s(szRealPath, L"\\\\.\\");
				wcscat_s(szRealPath, szVal);
			
				break;
			}
			//else
			//	msyslog(LOG_ERR, "VAR: failed to query aliased device port, error %u", lr);
		}

		SetupDiDestroyDeviceInfoList(hdi);
	}


	HANDLE  hCom = nullptr;
	hCom = ::CreateFile(szRealPath, GENERIC_READ | GENERIC_WRITE, FILE_SHARE_READ | FILE_SHARE_WRITE, NULL, OPEN_EXISTING, FILE_ATTRIBUTE_NORMAL | FILE_FLAG_OVERLAPPED, NULL );

	if (!hCom)
	{
		wprintf(L"There is no dev with %s", pszVID);
		return S_FALSE;
	}

	// we'll be signalled when 'carriage return' will be in buffer
	DCB dcb = { 0 };
    dcb.DCBlength = sizeof(dcb);
    assert( ::GetCommState( hCom, &dcb ) );
    

	dcb.fNull = TRUE;
	dcb.EvtChar = 0x0A;
	assert(::SetCommState( hCom, &dcb ) );

	//
	//-------------------------------------------------------------------------------------------
	COMMTIMEOUTS cto = { 0 };
	cto.ReadIntervalTimeout = MAXDWORD;		// this allows async ::ReadFile() to work immediately.
	assert(::SetCommTimeouts( hCom, &cto ) );

	//
	assert(::SetCommMask(hCom, EV_RXFLAG ) );

	CMHandle hEvt( ::CreateEvent( NULL, TRUE, FALSE, NULL ) );
	OVERLAPPED ovl = { 0 };
	ovl.hEvent = hEvt;
	bool fReenable = true;

	DWORD dwEventsAfterWait = 0;

	for( ; !g_fExit; /*::Sleep( 50 )*/ )
	{
		// слідкуємо асинхронно
		if( fReenable )
		{
			fReenable = false;

			auto Result = ::WaitCommEvent( hCom, &dwEventsAfterWait, &ovl );
			auto gle = ::GetLastError();
			assert( Result || gle == ERROR_IO_PENDING );
		}

		auto WaitResult = ::WaitForSingleObject( hEvt, INFINITE );	// ### add second event
		if( WaitResult == WAIT_TIMEOUT )
			wprintf_s( L"timed out" );
		else
		{
			assert( dwEventsAfterWait & EV_RXFLAG );

			OVERLAPPED ovlRead = { 0 };

			char pBuffer[4096] = { 0 };
			DWORD dwRead = 0;
			assert( ::ReadFile( hCom, pBuffer, sizeof(pBuffer) - 1, &dwRead, &ovlRead ) );
			DWORD gle;
			assert( (gle = ::GetLastError()) == ERROR_IO_PENDING || !gle );

			assert( !dwRead || pBuffer[dwRead-1] == 0x0A );

			if( dwRead )
			{
//				printf( "---------------------------------\r\n" );
//				printf( pBuffer );
			}
			else
				puts( "NOTHING!!!" );

			//
			// few lines possible got

			char seps[] = "\r\n";
			char *token = NULL;
			char *next_token = NULL;

			MyCoordinates Coords;
			// Establish string and get the first token:
			token = strtok_s( pBuffer, seps, &next_token );

			// While there are tokens in "string1" or "string2"
			while( token )
			{
				MiNMEA( token, Coords);

				token = strtok_s(NULL, seps, &next_token);
			}
			fReenable = true;

			if (Coords.m_fValid)
			{
				double dDeclination = 0;
				double dDeclinationError = 0;
				
				HRESULT hr;
				hr = Declination.GetDeclination(Coords.m_dLatitude, Coords.m_dLongitude, Coords.m_dHeight/1000, dDeclination, dDeclinationError);
				if (FAILED(hr))
					wprintf_s(L"Magnetic declination error");

				wprintf_s(L"Magnetic declination %f, Declination error %f \n", dDeclination, dDeclinationError);
			}
		}
	}

	return S_OK;
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

int wmain( int argc, LPWSTR argv[] )
{
	bool fPrnUsage = false;
	HRESULT hrMain = E_UNEXPECTED;

	WCHAR szErr[1024] = L"";

	std::wstring sAppName;

	do
	{
		::SetErrorMode( SEM_FAILCRITICALERRORS | SEM_NOOPENFILEERRORBOX );
		::SetConsoleCtrlHandler( ConsoleHandler, TRUE );

		hrMain = __GetModuleFileName( sAppName, NULL, NULL, MGMN_ONLYNAME_ );
		if ( FAILED( hrMain ) )
			break;

		if ( argc != 1 )
		{
			hrMain = E_INVALIDARG;
			fPrnUsage = true;
			break;
		}

		hrMain = DoMain();
		if ( FAILED( hrMain ) )
			break;
	}
	while ( 0 );

	if ( fPrnUsage )
	{
		LPCWSTR pszApp = (WCHAR*)sAppName.c_str();
		if ( !pszApp[0] )
			pszApp = ::PathFindFileName( argv[0] );
		
		if ( !pszApp || !pszApp[0] )
			pszApp = L"exe";

		wprintf_s( L"usage: %s, without any cmdline arguments\n", pszApp );
	}
	else if ( FAILED( hrMain ) )
	{
		wprintf_s( L"finished with error: %d\n", szErr );
	}

	return hrMain;
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
