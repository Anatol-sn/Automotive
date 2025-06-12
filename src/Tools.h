#ifndef TOOLS_H
#define TOOLS_H

#pragma once

#include <winnt.h>
#include <Psapi.h>

#include <shlwapi.h>
#pragma comment( lib, "Shlwapi.lib")

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

typedef enum MGETMODULEFILENAMEFLAGS_
{
	MGMN_FULL_ = 0x0,				// default
	MGMN_TRIMEXT_ = 0x1,			// без розширення і попередньої точки (якщо таке є)
	MGMN_TRIMNAME_ = 0x2,			// без імені файлу (і розширення), останній символ - '\'
	MGMN_ONLYNAME_ = 0x3,			// тільки ім'я файлу (і розширення)
	MGMN_ONLYNAME_TRIMEXT_ = 0x4,	// тільки ім'я файлу, без розширення
} MGETMODULEFILENAMEFLAGS_;
//---------------------------------------------------------------------------

inline HRESULT MGetLastError()
{
	const DWORD dwErr = ::GetLastError();

	if (ERROR_SUCCESS == dwErr)
	{
		assert(!L"invalid condition");
		return E_FAIL;
	}

	return HRESULT_FROM_WIN32(dwErr);
}
//---------------------------------------------------------------------------

#define M_DECLARE_NOCOPY( cls )				\
	cls( const cls& ) = delete;				\
	cls& operator=( const cls& ) = delete;

#define M_DECLARE_NOMOVE( cls )				\
	cls( cls&& ) = delete;					\
	cls& operator=( cls&& ) = delete;

#define M_DECLARE_DEFCOPY( cls )			\
	cls( const cls& ) = default;			\
	cls& operator=( const cls& ) = default;

#define M_DECLARE_DEFMOVE( cls )			\
	cls( cls&& ) = default;					\
	cls& operator=( cls&& ) = default;

//---------------------------------------------------------------------------
#define mdecl_noexcept noexcept
//---------------------------------------------------------------------------
struct m_declare_move
{};
#define m_confirm_move m_declare_move()
//---------------------------------------------------------------------------
#define mmove std::move

class CMHandle
{
	M_DECLARE_NOCOPY(CMHandle)

protected:
	HANDLE m_h;

public:
	CMHandle()
		mdecl_noexcept
		: m_h(nullptr)
	{}

	~CMHandle()
	{
		zClose(m_h);
		//m_h = nullptr;
	}

	explicit CMHandle(HANDLE h)
		mdecl_noexcept
		: m_h(h)
	{}

	CMHandle(CMHandle& r, m_declare_move)
		mdecl_noexcept
		: m_h(r.Detach())
	{}

	CMHandle(CMHandle&& r)
		mdecl_noexcept
		: m_h(r.Detach())
	{}

	CMHandle& operator=(CMHandle&& r)
	{
		TakeAway(mmove(r));
		return (*this);
	}

	operator HANDLE() const mdecl_noexcept
	{
		return m_h;
	}

	HANDLE Get() const mdecl_noexcept
	{
		return m_h;
	}

	const HANDLE* GetAddr() const mdecl_noexcept { return &m_h; }
	HANDLE* GetAddr() mdecl_noexcept { return &m_h; }

	bool IsValid() const mdecl_noexcept
	{
		return !!m_h;
	}

	void Attach(HANDLE h) mdecl_noexcept
	{
		assert(!m_h);
		m_h = h;
	}

	HANDLE Detach() mdecl_noexcept
	{
		HANDLE h = m_h;
		m_h = nullptr;
		return h;
	}

	void Close()
	{
		HANDLE h = m_h;
		m_h = nullptr;
		zClose(h);
	}

	void Swap(CMHandle& r) mdecl_noexcept
	{
		std::swap(m_h, r.m_h);
	}

	void TakeAway(CMHandle&& r)
	{
		// take into account arg == this
		HANDLE hNew = r.Detach();
		HANDLE hOld = m_h;
		m_h = hNew;
		zClose(hOld);
	}

protected:
	static inline void zClose(HANDLE h)
	{
		if (h)
			assert(::CloseHandle(h));
	}
};
//---------------------------------------------------------------------------

HRESULT __GetModuleFileName(std::wstring& rsOut, HMODULE hModule, HANDLE hProcess, DWORD dwFlags)
{
	WCHAR pszBuf[MAX_PATH];

	while (1)
	{
		DWORD dwLen;
		if (hProcess)
			dwLen = ::GetModuleFileNameEx(hProcess, hModule, pszBuf, (DWORD)MAX_PATH);
		else
			dwLen = ::GetModuleFileName(hModule, pszBuf, (DWORD)MAX_PATH);
		if (!dwLen)
			return MGetLastError();							\

			if (MAX_PATH > dwLen)
				break;
	}

	LPWSTR pszResult = pszBuf;

	if (MGMN_FULL_ == dwFlags)
		;
	else if (MGMN_TRIMEXT_ == dwFlags)
	{
		LPWSTR psz = ::PathFindExtension(pszBuf);
		if (!psz)
			return E_FAIL;

		if (L'.' == *psz)
			*psz = 0;
		else
			if (0 != *psz)
				return E_FAIL;
	}
	else if (MGMN_TRIMNAME_ == dwFlags)
	{
		LPWSTR psz = ::PathFindFileName(pszBuf);
		if (!psz)
			return E_FAIL;
		*psz = 0;
	}
	else if (MGMN_ONLYNAME_ == dwFlags)
	{
		pszResult = ::PathFindFileName(pszBuf);
		if (!pszResult)
			return E_FAIL;
	}
	else if (MGMN_ONLYNAME_TRIMEXT_ == dwFlags)
	{
		pszResult = ::PathFindFileName(pszBuf);
		if (!pszResult)
			return E_FAIL;

		LPWSTR psz = ::PathFindExtension(pszResult);
		if (!psz)
			return E_FAIL;

		if (L'.' == *psz)
			*psz = 0;
		else
			if (0 != *psz)
				return E_FAIL;
	}
	else
		return E_INVALIDARG;

	rsOut = pszResult;
	return S_OK;
}
//---------------------------------------------------------------------------

#endif //TOOLS_H