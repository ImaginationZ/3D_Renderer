// CS580HWView.cpp : implementation of the CCS580HWView class
//

#include "stdafx.h"
#include "CS580HW.h"

#include "CS580HWDoc.h"
#include "CS580HWView.h"
#include "RotateDlg.h"
#include "TranslateDlg.h"
#include "ScaleDlg.h"

#include "disp.h"
#include "ApplicationFinal.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

/////////////////////////////////////////////////////////////////////////////
// CCS580HWView

IMPLEMENT_DYNCREATE(CCS580HWView, CView)

BEGIN_MESSAGE_MAP(CCS580HWView, CView)
	//{{AFX_MSG_MAP(CCS580HWView)
	ON_COMMAND(IDM_RENDER, OnRender)
	ON_COMMAND(IDM_ROTATE, OnRotate)
	ON_COMMAND(IDM_TRANSLATE, OnTranslate)
	ON_COMMAND(IDM_SCALE, OnScale)
	ON_COMMAND(ID_MATERIAL_REFLECTIVE, OnMaterialReflective)
	//}}AFX_MSG_MAP
	ON_UPDATE_COMMAND_UI(ID_MATERIAL_REFLECTIVE, &CCS580HWView::OnUpdateMaterialReflective)
	ON_COMMAND(ID_MATERIAL_GLASS, &CCS580HWView::OnMaterialGlass)
	ON_UPDATE_COMMAND_UI(ID_MATERIAL_GLASS, &CCS580HWView::OnUpdateMaterialGlass)
	ON_COMMAND(ID_MATERIAL_WATER, &CCS580HWView::OnMaterialWater)
	ON_UPDATE_COMMAND_UI(ID_MATERIAL_WATER, &CCS580HWView::OnUpdateMaterialWater)
	ON_COMMAND(ID_MATERIAL_DIAMOND, &CCS580HWView::OnMaterialDiamond)
	ON_UPDATE_COMMAND_UI(ID_MATERIAL_DIAMOND, &CCS580HWView::OnUpdateMaterialDiamond)
	ON_COMMAND(ID_MATERIAL_AIR, &CCS580HWView::OnMaterialAir)
	ON_UPDATE_COMMAND_UI(ID_MATERIAL_AIR, &CCS580HWView::OnUpdateMaterialAir)
END_MESSAGE_MAP()

/////////////////////////////////////////////////////////////////////////////
// CCS580HWView construction/destruction

CCS580HWView::CCS580HWView()
{
	// TODO: add construction code here
	m_pApplication = NULL;
	m_isReflective = false;
	m_isWater = false;
	m_isGlass = false;
	m_isDiamond = false;
	m_isAir = false;
}

CCS580HWView::~CCS580HWView()
{
	if(m_pApplication != NULL)
	{
		delete m_pApplication;
	}
}

BOOL CCS580HWView::PreCreateWindow(CREATESTRUCT& cs)
{
	// TODO: Modify the Window class or styles here by modifying
	//  the CREATESTRUCT cs

	return CView::PreCreateWindow(cs);
}

/////////////////////////////////////////////////////////////////////////////
// CCS580HWView drawing

void CCS580HWView::OnDraw(CDC* pDC)
{
	CCS580HWDoc* pDoc = GetDocument();
	ASSERT_VALID(pDoc);
	// TODO: add draw code for native data here
	if(m_pApplication != NULL)
		DrawFrameBuffer(pDC);
}

/////////////////////////////////////////////////////////////////////////////
// CCS580HWView diagnostics

#ifdef _DEBUG
void CCS580HWView::AssertValid() const
{
	CView::AssertValid();
}

void CCS580HWView::Dump(CDumpContext& dc) const
{
	CView::Dump(dc);
}

CCS580HWDoc* CCS580HWView::GetDocument() // non-debug version is inline
{
	ASSERT(m_pDocument->IsKindOf(RUNTIME_CLASS(CCS580HWDoc)));
	return (CCS580HWDoc*)m_pDocument;
}
#endif //_DEBUG

/////////////////////////////////////////////////////////////////////////////
// CCS580HWView message handlers

void CCS580HWView::OnRender() 
{
	// TODO: Add your command handler code here

	// Call renderer 

	// Application Final
	if(m_pApplication != NULL)
		((ApplicationFinal *)m_pApplication)->Render();
	else 
		AfxMessageBox("Application was not allocated\n");

	// Set window size
	CRect clientRect, windowRect;
	int x_offset, y_offset;

	GetClientRect(&clientRect);
	AfxGetMainWnd()->GetWindowRect(&windowRect);
	
	x_offset = windowRect.Width() - clientRect.Width();
	y_offset = windowRect.Height() - clientRect.Height();

	AfxGetMainWnd()->SetWindowPos(NULL, 0, 0, x_offset+m_pApplication->m_nWidth, y_offset+m_pApplication->m_nHeight, NULL/*,SWP_SHOWWINDOW*/);

	Invalidate(true);	
}

void CCS580HWView::DrawFrameBuffer(CDC *pDC)
{
	if(m_pApplication->m_pFrameBuffer == NULL)
    {
        return;
    }

    HDC hdc;
    hdc = ::CreateCompatibleDC(pDC->m_hDC);
	HBITMAP m_bitmap;

    // Display the current image
    char buffer[sizeof(BITMAPINFO)];
    BITMAPINFO* binfo = (BITMAPINFO*)buffer;
    memset(binfo, 0, sizeof(BITMAPINFOHEADER));
    binfo->bmiHeader.biSize = sizeof(BITMAPINFOHEADER);
    
	// Create the bitmap
    BITMAPINFOHEADER* bih = &binfo->bmiHeader;
	bih->biBitCount = 3*8;
	bih->biWidth = m_pApplication->m_nWidth;
	bih->biHeight = m_pApplication->m_nHeight;
    bih->biPlanes = 1;
    bih->biCompression = BI_RGB;
    bih->biSizeImage = 0;
    
    m_bitmap = CreateDIBSection(hdc, binfo, 0, 0, 0, DIB_RGB_COLORS);

    int colors = DIB_RGB_COLORS;
    
    ::SelectObject(hdc, m_bitmap);
	binfo->bmiHeader.biBitCount = 0;
	GetDIBits(hdc, m_bitmap, 0, 0, 0, binfo, colors);
    binfo->bmiHeader.biBitCount = 24;
    binfo->bmiHeader.biHeight = -abs(binfo->bmiHeader.biHeight);
    SetDIBits(hdc, m_bitmap, 0, m_pApplication->m_nHeight, m_pApplication->m_pFrameBuffer, binfo, colors);

    ::SetStretchBltMode(pDC->m_hDC, COLORONCOLOR);
    CRect client;
    GetClientRect(&client);
    ::BitBlt(pDC->m_hDC, 0, 0, m_pApplication->m_nWidth, m_pApplication->m_nHeight, 
        hdc, 0, 0, SRCCOPY);
    ::DeleteDC(hdc);
    DeleteObject(m_bitmap); 
}

void CCS580HWView::OnInitialUpdate() 
{
	CView::OnInitialUpdate();
	
	// TODO: Add your specialized code here and/or call the base class

	// Assign Application Final
	if(m_pApplication == NULL)
	{
		m_pApplication = new ApplicationFinal;
	}
	
	// Initialize and begin renderer
	((ApplicationFinal *)m_pApplication)->Initialize();
}

// Callback function for rotation  
void CCS580HWView::OnRotate() 
{
	// TODO: Add your command handler code here
	CRotateDlg dlg;
	GzInput* input;

	if(m_pApplication == NULL) return;

	input = m_pApplication->m_pUserInput;
	if(input == NULL) return;

	// Initialize
	input->rotation[0] = input->rotation[1] = input->rotation[2] = 0;
	dlg.Initialize(input->rotation[0], input->rotation[1], input->rotation[2]);

	if(dlg.DoModal() == IDOK)
	{
		// Update input rotation value
		input->rotation[dlg.m_nAxis] = dlg.m_fRot;

		// Accumulate matrix
		for (int aaPass = 0; aaPass < AAKERNEL_SIZE; aaPass++) {
			GzMatrix	rotMat = 
			{ 
				1.0,	0.0,	0.0,	0.0, 
				0.0,	1.0,	0.0,	0.0, 
				0.0,	0.0,	1.0,	0.0, 
				0.0,	0.0,	0.0,	1.0 
			};
			switch(dlg.m_nAxis)
			{
			case 0 :
				// Create matrix for Rot X
				GzRotXMat(input->rotation[0], rotMat);
				break;
			case 1:
				// Create matrix for Rot Y
				GzRotYMat(input->rotation[1], rotMat);
				break;
			case 2:
				// Create matrix for Rot Z
				GzRotZMat(input->rotation[2], rotMat);
				break;
			}
			GzPushMatrix(m_pApplication->m_pAARenders[aaPass], rotMat); 
		}
	}
}

// Callback function for Translation
void CCS580HWView::OnTranslate() 
{
	// TODO: Add your command handler code here
	CTranslateDlg dlg;
	GzInput* input;

	if(m_pApplication == NULL) return;

	input = m_pApplication->m_pUserInput;
	if(input == NULL) return;

	// Initialize
	input->translation[0] = input->translation[1] = input->translation[2] = 0;
	dlg.Initialize(input->translation[0], input->translation[1], input->translation[2]);

	if(dlg.DoModal() == IDOK)
	{
		// Update input translation value
		input->translation[0] = dlg.m_fTx; input->translation[1] = dlg.m_fTy; input->translation[2] = dlg.m_fTz;

		// Accumulate matrix
		for (int aaPass = 0; aaPass < AAKERNEL_SIZE; aaPass++) {
			GzMatrix trxMat = 
			{ 
				1.0,	0.0,	0.0,	0.0, 
				0.0,	1.0,	0.0,	0.0, 
				0.0,	0.0,	1.0,	0.0, 
				0.0,	0.0,	0.0,	1.0 
			};
			GzTrxMat(input->translation, trxMat);
			GzPushMatrix(m_pApplication->m_pAARenders[aaPass], trxMat); 
		}
	}
}

// Callback function for Scaling
void CCS580HWView::OnScale() 
{
	// TODO: Add your command handler code here
	CScaleDlg dlg;
	GzInput* input;

	if(m_pApplication == NULL) return;

	input = m_pApplication->m_pUserInput;
	if(input == NULL) return;

	// Initialize
	input->scale[0] = input->scale[1] = input->scale[2] = 1;
	dlg.Initialize(input->scale[0], input->scale[1], input->scale[2]);

	if(dlg.DoModal() == IDOK)
	{
		// Update input scale value
		input->scale[0] = dlg.m_fSx; input->scale[1] = dlg.m_fSy; input->scale[2] = dlg.m_fSz;

		// Accumulate matrix
		for (int aaPass = 0; aaPass < AAKERNEL_SIZE; aaPass++) {
			GzMatrix scaleMat = 
			{ 
				1.0,	0.0,	0.0,	0.0, 
				0.0,	1.0,	0.0,	0.0, 
				0.0,	0.0,	1.0,	0.0, 
				0.0,	0.0,	0.0,	1.0 
			};
			GzScaleMat(input->scale, scaleMat);
			GzPushMatrix(m_pApplication->m_pAARenders[aaPass], scaleMat); 
		}
	}
}


void CCS580HWView::OnMaterialReflective()
{
	bool temp = m_isReflective;
	ResetCheckboxes();

	m_isReflective = !temp;
	m_pApplication->SetReflective(m_isReflective);
}


void CCS580HWView::OnUpdateMaterialReflective(CCmdUI *pCmdUI)
{
	pCmdUI->SetCheck(m_isReflective);
}


void CCS580HWView::OnMaterialGlass()
{
	bool temp = m_isGlass;
	ResetCheckboxes();

	m_isGlass = !temp;
	m_pApplication->SetRefractive(m_isGlass);
	m_pApplication->SetRefractiveIndex(1.5);
}


void CCS580HWView::OnUpdateMaterialGlass(CCmdUI *pCmdUI)
{
	pCmdUI->SetCheck(m_isGlass);
}


void CCS580HWView::OnMaterialWater()
{
	bool temp = m_isWater;
	ResetCheckboxes();

	m_isWater = !temp;
	m_pApplication->SetRefractive(m_isWater);
	m_pApplication->SetRefractiveIndex(1.325);
}


void CCS580HWView::OnUpdateMaterialWater(CCmdUI *pCmdUI)
{
	pCmdUI->SetCheck(m_isWater);
}


void CCS580HWView::OnMaterialDiamond()
{
	bool temp = m_isDiamond;
	ResetCheckboxes();

	m_isDiamond = !temp;
	m_pApplication->SetRefractive(m_isDiamond);
	m_pApplication->SetRefractiveIndex(2.418);
}


void CCS580HWView::OnUpdateMaterialDiamond(CCmdUI *pCmdUI)
{
	pCmdUI->SetCheck(m_isDiamond);
}


void CCS580HWView::OnMaterialAir()
{
	bool temp = m_isAir;
	ResetCheckboxes();

	m_isAir = !temp;
	m_pApplication->SetRefractive(m_isAir);
	m_pApplication->SetRefractiveIndex(1.0);
}


void CCS580HWView::OnUpdateMaterialAir(CCmdUI *pCmdUI)
{
	pCmdUI->SetCheck(m_isAir);
}

void CCS580HWView::ResetCheckboxes() 
{
	m_isReflective = false;
	m_isWater = false;
	m_isGlass = false;
	m_isDiamond = false;
	m_isAir = false;

	m_pApplication->SetReflective(false);
	m_pApplication->SetRefractive(false);
	m_pApplication->SetRefractiveIndex(0);
}
