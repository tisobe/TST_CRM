; plot Chandra orbital proton fluence

; Robert Cameron
; September 2000
; David Morris
; December 2000 : updated to new CRM_p_mn.dat format and year rollover
; RAC, Feb 2001: updated to new CRM data file format.

crmdir = '/data/mta4/proj/rac/ops/CRM/'

; set up colour vectors

red = bytarr(256)+255
green = red
blue = red
  red[0:6]=[0, 255, 0,   0,   0,   255, 255]
green[0:6]=[0, 0,   255, 0,   255, 255, 0]
 blue[0:6]=[0, 0,   0,   255, 255, 0,   255]
TVLCT, red,green,blue

; create a 100-year vector of day offsets from 1998

ylen = [366,365,365,365]
d98 = lonarr(100)
for i = 1,99 do d98(i) = d98(i-1) + ylen((1998 + i-1) mod 4)

; get the current time, in seconds since 1998

t1998 = systime(1) - 883612800L;
timestring = "   (Plot generation time = " + systime() + ")";

; get the current CRM summary data

spawn,"perl -ane '@t = split; print ""$t[-1]\n""' /data/mta4/proj/rac/ops/CRM/CRMsummary.dat",sum
kp = sum[1]
ace = double(sum[2])
fluence = double(sum[12])
afluence = double(sum[13])

; convert and read the FPSI history file

f1 = crmdir+"FPHIST-2001.dat"
sif = crmdir+"fphist.dat"
spawn,"perl -pe 's/:/ /g; s/H/1/; s/-S/1/; s/-I/2/; tr/A-Z//d' "+f1+" >"+sif
rac_arread,sif,ysi,dsi,hsi,msi,ssi,si,obssi
dsi = d98(ysi-1998) + dsi + hsi/24. + msi/1400. + ssi/86400.

nsi = n_elements(dsi)
dsi1 = shift(dsi,-1)
dsi1(nsi-1) = 99999.0
ok = where(si gt 2,nok)
if nok gt 0 then si(ok) = si(ok)-8

; convert and read the OTG history file

f1 = crmdir+"GRATHIST-2001.dat"
otgf = crmdir+"grathist.dat"
spawn,"perl -pe 's/:/ /g; s/I/1/g; s/O/0/g; tr/A-Z-//d' "+f1+" >"+otgf
rac_arread,otgf,yotg,dotg,hotg,motg,sotg,hetg,letg,obsotg
dotg = d98(yotg-1998) + dotg + hotg/24. + motg/1400. + sotg/86400.

notg = n_elements(dotg)
dotg1 = shift(dotg,-1)
dotg1(notg-1) = 99999.0
hidx=where(hetg gt 0,nhetg)
if nhetg gt 0 then begin
  dhetg = dotg(hidx)
  dhetg1 = dotg1(hidx)
endif
lidx=where(letg gt 0,nletg)
if nletg gt 0 then begin
  dletg = dotg(lidx)
  dletg1 = dotg1(lidx)
endif

; read the DSN contact schedule

filename = '/data/mta4/proj/rac/ops/ephem/DSN.sch'
spawn,'wc -l '+filename,info
info=str_sep(strcompress(info(0)),' ')
nlines=long(info(1))
r = {byr:0,bot:0.0,eyr:0,eot:0.0,txt:''}
r=replicate(r,nlines)
openr,nlun,filename,/get_lun
readf,nlun,r
bot = d98(r.byr-1998) + r.bot
eot = d98(r.eyr-1998) + r.eot
mot = (bot + eot)/2

; read the GSM, GSE ephemerides

rac_arread,'/data/mta4/proj/rac/ops/ephem/PE.EPH.gsme',t,r,mlat,mlon,elat,elon,fy,mon,day,hour,min,sec
r = r/1e3

; read the CRM regions and fluxes

farr = findfile(crmdir+'CRM_p.dat*0',count=nfarr)
kpext = string(nint(kp*10.0),format='(i2.2)')
farr = [crmdir+'CRM_p.dat'+kpext,farr]
nfarr = nfarr + 1
k = dblarr(n_elements(t),nfarr)
reg = intarr(n_elements(t),nfarr)
for i = 0,nfarr-1 do begin
  rac_arread,farr(i),t,r0,fmn,f95,f50,fsd
  reg(*,i) = r0
  k(*,i) = fmn
endfor

; make time arrays

timeconv,t,y,d,h,m,s,y0=1998
y = fix(y)-1998
d = d + d98(y)
de = rebin(d,n_elements(d),nfarr)

; truncate the arrays to cover only one orbit

past = where(t lt t1998,npast);
if npast gt 0 then begin
  k(past,*) = 0
  idx1 = (npast - 10) > 0
  idx2 = (idx1 + 8000) < (n_elements(t)-1)
  d=d[idx1:idx2]
  y=y[idx1:idx2]
  fy = fy[idx1:idx2]
  r=r[idx1:idx2]
  k=k[idx1:idx2,*]
  reg=reg[idx1:idx2,*]
  de=de[idx1:idx2,*]
endif

; find the perigees

delr1 = r - shift(r,-1)
delr2 = shift(delr1,-1)
pg = intarr(n_elements(r))+1
ok = where(delr1 ge 0 and delr2 le 0,nok)
if nok gt 0 then pg(ok) = 0

; truncate the SI and OTG history arrays to the time range

nd = n_elements(d)
if nsi gt 0 then $
ok = where(dsi le d(nd-1) and dsi1 ge d(0),nsi)
if nsi gt 0 then begin
  dsi = dsi(ok)
  dsi1 = dsi1(ok)
  si = si(ok)
endif
if nletg gt 0 then $
ok = where(dletg le d(nd-1) and dletg1 ge d(0),nletg)
if nletg gt 0 then begin
  dletg = dletg(ok)
  dletg1 = dletg1(ok)
endif
if nhetg gt 0 then $
ok = where(dhetg le d(nd-1) and dhetg1 ge d(0),nhetg)
if nhetg gt 0 then begin
  dhetg = dhetg(ok)
  dhetg1 = dhetg1(ok)
endif

; create FPSI and OTG attenuation factor vectors

siv = fltarr(nd)
siaf = siv+1
letgaf = siaf
hetgaf = siaf
for i=0,nd-1 do begin
  di = d(i)
  for j=0,nletg-1 do if (di gt dletg(j) and di le dletg1(j)) then letgaf(i) = 0.5
  for j=0,nhetg-1 do if (di gt dhetg(j) and di le dhetg1(j)) then hetgaf(i) = 0.2
  for j=0,nsi-1 do if (di gt dsi(j) and di le dsi1(j)) then siv(i) = si(j)
endfor
ok=where(siv gt 2, nok)
if nok gt 0 then siaf(ok) = 0.0
af = siaf*letgaf*hetgaf

; set the current fluence at the current time in the k array
; and integrate either the CRM or current ACE fluence

reg1= where(reg eq 1,nreg1)
reg2= where(reg eq 2,nreg2)
reg3= where(reg eq 3,nreg3)
kreg1= where(reg[*,0] eq 1,nkreg1)
kreg2= where(reg[*,0] eq 2,nkreg2)
kreg3= where(reg[*,0] eq 3,nkreg3)

if nreg1 gt 0 then k(reg1)=ace
if nreg2 gt 0 then k(reg2)=k(reg2)+2*ace
k=k*300
ki = k
ka = k
kia = ka
ka = ka * rebin(af,nd,nfarr)
ki(0,*) = fluence
kia(0,*) = afluence
for i = 1,n_elements(r)-1 do begin  
  ki(i,*) = ki(i-1,*)*pg(i) + k(i,*)
  kia(i,*) = kia(i-1,*)*pg(i) + ka(i,*)
endfor

; make the plots

start_JD = julday(1,0,1998)
dummy = LABEL_DATE(DATE_FORMAT='%M-%D',offset=start_JD)

!x.style=1
!p.charsize=2
!y.omargin=[2,0]

; the EXTERNAL fluence plot

!p.multi=[0,1,5]

; plot altitude and DSN contact schedule

plot, d, r, ytit='Altitude (Mm)',ymargin=[0,2],xmargin=[10,4],yrange=[0,150],$
  tit='EXTERNAL CRM Proton Fluence, with DSN and SI schedules',yticks=3,$
  xminor=12,xtickn=[' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ']
for i=0,n_elements(bot)-1 do oplot, [bot(i),eot(i)],[75,75],thick=70
if nkreg1 gt 0 then oplot, d(kreg1),r(kreg1),psym=3,color=4
if nkreg2 gt 0 then oplot, d(kreg2),r(kreg2),psym=3,color=6
if nkreg3 gt 0 then oplot, d(kreg3),r(kreg3),psym=3,color=5

; plot SI and OTG configs
; key: ACIS-S = 1, ACIS-I = 2, HRC-S = 3, HRC-I = 4

plot, [d(0),d(n_elements(d)-1)]-d98(y(0)),[0,0],yra=[0,7],ymar=[0,0],xmar=[10,4],xminor=12,$
  yticks=7,yminor=-1, ytickname=[' ','ACIS-S','ACIS-I','HRC-S','HRC-I','HETG','LETG',' ']
for i=0,n_elements(mot)-1 do oplot, [mot(i),mot(i)]-d98(y(0)),[0,7]
for i=1,6 do oplot, [d(0),d(n_elements(d)-1)]-d98(y(0)),[i,i]
for i=0,nsi-1 do oplot,[dsi(i),dsi1(i)]-d98(y(0)),[si(i),si(i)],thick=4,color=1
if nhetg gt 0 then for i=0,nhetg-1 do oplot,[dhetg(i),dhetg1(i)]-d98(y(0)),[5,5],thick=4,color=2
if nletg gt 0 then for i=0,nletg-1 do oplot,[dletg(i),dletg1(i)]-d98(y(0)),[6,6],thick=4,color=2

; plot CRM model

plot_io, d,ki(*,0), yra=[min(ki)>1e7,max(ki)<1e12>1e10],XTICKFORMAT='LABEL_DATE',xminor=12,$
xtitle='UTC Date, Day of Year' + timestring,xmargin=[10,4],ymargin=[-14.5,0],ytit='Proton Fluence (p/cm2-sr-MeV)'
if nreg1 gt 0 then oplot, de(reg1),ki(reg1),psym=3,color=4
if nreg2 gt 0 then oplot, de(reg2),ki(reg2),psym=3,color=6
if nreg3 gt 0 then oplot, de(reg3),ki(reg3),psym=3,color=5
oplot, d,ki(*,0),color=2
for i=0,n_elements(mot)-1 do oplot, [mot(i),mot(i)],10^!y.crange
oplot, [d(0),d(n_elements(d)-1)],[3e9,3e9],color=1

out=tvrd()
WRITE_GIF, crmdir+'crmpl.gif', out,red,green,blue

; the ATTENUATED fluence plot

!p.multi=[0,1,5]

; plot altitude and DSN contact schedule

plot, d, r, ytit='Altitude (Mm)',ymargin=[0,2],xmargin=[10,4],yrange=[0,150],$
  tit='ATTENUATED CRM Proton Fluence, with DSN and SI schedules',yticks=3,$
  xminor=12,xtickn=[' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ']
for i=0,n_elements(bot)-1 do oplot, [bot(i),eot(i)],[75,75],thick=70
if nkreg1 gt 0 then oplot, d(kreg1),r(kreg1),psym=3,color=4
if nkreg2 gt 0 then oplot, d(kreg2),r(kreg2),psym=3,color=6
if nkreg3 gt 0 then oplot, d(kreg3),r(kreg3),psym=3,color=5

; plot SI and OTG configs
; key: ACIS-S = 1, ACIS-I = 2, HRC-S = 3, HRC-I = 4

plot, [d(0),d(n_elements(d)-1)]-d98(y(0)),[0,0],yra=[0,7],ymar=[0,0],xmar=[10,4],xminor=12,$
  yticks=7,yminor=-1, ytickname=[' ','ACIS-S','ACIS-I','HRC-S','HRC-I','HETG','LETG',' ']
for i=0,n_elements(mot)-1 do oplot, [mot(i),mot(i)]-d98(y(0)),[0,7]
for i=1,6 do oplot, [d(0),d(n_elements(d)-1)]-d98(y(0)),[i,i]
for i=0,nsi-1 do oplot,[dsi(i),dsi1(i)]-d98(y(0)),[si(i),si(i)],thick=4,color=1
if nhetg gt 0 then for i=0,nhetg-1 do oplot,[dhetg(i),dhetg1(i)]-d98(y(0)),[5,5],thick=4,color=2
if nletg gt 0 then for i=0,nletg-1 do oplot,[dletg(i),dletg1(i)]-d98(y(0)),[6,6],thick=4,color=2

; plot CRM model

plot_io, d,kia(*,0), yra=[min(kia)>1e7,max(kia)<1e12>1e10],XTICKFORMAT='LABEL_DATE',xminor=12,$
xtitle='UTC Date, Day of Year' + timestring,xmargin=[10,4],ymargin=[-14.5,0],ytit='Proton Fluence (p/cm2-sr-MeV)'
if nreg1 gt 0 then oplot, de(reg1),kia(reg1),psym=3,color=4
if nreg2 gt 0 then oplot, de(reg2),kia(reg2),psym=3,color=6
if nreg3 gt 0 then oplot, de(reg3),kia(reg3),psym=3,color=5
oplot, d,kia(*,0),color=2
for i=0,n_elements(mot)-1 do oplot, [mot(i),mot(i)],10^!y.crange
oplot, [d(0),d(n_elements(d)-1)],[3e9,3e9],color=1

out=tvrd()
WRITE_GIF, crmdir+'crmplatt.gif', out,red,green,blue

; the EXTERNAL flux plot
;
;!p.multi=[0,1,5]
;
; plot altitude and DSN contact schedule
;
;plot, d, r, ytit='Altitude (Mm)',ymargin=[0,2],xmargin=[10,4],yrange=[0,150],$
;  tit='Future EXTERNAL Proton Flux, with DSN and SI schedules',yticks=3,$
;  xminor=12,xtickn=[' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ']
;for i=0,n_elements(bot)-1 do oplot, [bot(i),eot(i)],[75,75],thick=70
;if nkreg1 gt 0 then oplot, d(kreg1),r(kreg1),psym=3,color=4
;if nkreg2 gt 0 then oplot, d(kreg2),r(kreg2),psym=3,color=6
;if nkreg3 gt 0 then oplot, d(kreg3),r(kreg3),psym=3,color=5
;
; plot SI and OTG configs
; key: ACIS-S = 1, ACIS-I = 2, HRC-S = 3, HRC-I = 4
;
;plot, [d(0),d(n_elements(d)-1)]-d98(y(0)),[0,0],yra=[0,7],ymar=[0,0],xmar=[10,4],xminor=12,$
;  yticks=7,yminor=-1, ytickname=[' ','ACIS-S','ACIS-I','HRC-S','HRC-I','HETG','LETG',' ']
;for i=0,n_elements(mot)-1 do oplot, [mot(i),mot(i)]-d98(y(0)),[0,7]
;for i=1,6 do oplot, [d(0),d(n_elements(d)-1)]-d98(y(0)),[i,i]
;for i=0,nsi-1 do oplot,[dsi(i),dsi1(i)]-d98(y(0)),[si(i),si(i)],thick=4,color=1
;if nhetg gt 0 then for i=0,nhetg-1 do oplot,[dhetg(i),dhetg1(i)]-d98(y(0)),[5,5],thick=4,color=2
;if nletg gt 0 then for i=0,nletg-1 do oplot,[dletg(i),dletg1(i)]-d98(y(0)),[6,6],thick=4,color=2
;
; plot CRM model
;
;plot_io, d,k(*,0), yra=[min(k)>1,max(k)<1e10],XTICKFORMAT='LABEL_DATE',xminor=12,$
;xtitle='UTC Date, Day of Year',xmargin=[10,4],ymargin=[-14.5,0],ytit='Proton Flux (p/s-cm2-sr-MeV)'
;if nreg1 gt 0 then oplot, de(reg1),k(reg1),psym=3,color=4
;if nreg2 gt 0 then oplot, de(reg2),k(reg2),psym=3,color=6
;if nreg3 gt 0 then oplot, de(reg3),k(reg3),psym=3,color=5
;oplot, d,k(*,0),color=2
;for i=0,n_elements(mot)-1 do oplot, [mot(i),mot(i)],10^!y.crange
;oplot, [d(0),d(n_elements(d)-1)],[5e4,5e4],color=1
;
;out=tvrd()
;WRITE_GIF, crmdir+'crmplx.gif', out,red,green,blue
;
; the ATTENUATED flux plot
;
;!p.multi=[0,1,5]
;
; plot altitude and DSN contact schedule
;
;plot, d, r, ytit='Altitude (Mm)',ymargin=[0,2],xmargin=[10,4],yrange=[0,150],$
;  tit='Future ATTENUATED Proton Flux, with DSN and SI schedules',yticks=3,$
;  xminor=12,xtickn=[' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ']
;for i=0,n_elements(bot)-1 do oplot, [bot(i),eot(i)],[75,75],thick=70
;if nkreg1 gt 0 then oplot, d(kreg1),r(kreg1),psym=3,color=4
;if nkreg2 gt 0 then oplot, d(kreg2),r(kreg2),psym=3,color=6
;if nkreg3 gt 0 then oplot, d(kreg3),r(kreg3),psym=3,color=5
;
; plot SI and OTG configs
; key: ACIS-S = 1, ACIS-I = 2, HRC-S = 3, HRC-I = 4
;
;plot, [d(0),d(n_elements(d)-1)]-d98(y(0)),[0,0],yra=[0,7],ymar=[0,0],xmar=[10,4],xminor=12,$
;  yticks=7,yminor=-1, ytickname=[' ','ACIS-S','ACIS-I','HRC-S','HRC-I','HETG','LETG',' ']
;for i=0,n_elements(mot)-1 do oplot, [mot(i),mot(i)]-d98(y(0)),[0,7]
;for i=1,6 do oplot, [d(0),d(n_elements(d)-1)]-d98(y(0)),[i,i]
;for i=0,nsi-1 do oplot,[dsi(i),dsi1(i)]-d98(y(0)),[si(i),si(i)],thick=4,color=1
;if nhetg gt 0 then for i=0,nhetg-1 do oplot,[dhetg(i),dhetg1(i)]-d98(y(0)),[5,5],thick=4,color=2
;if nletg gt 0 then for i=0,nletg-1 do oplot,[dletg(i),dletg1(i)]-d98(y(0)),[6,6],thick=4,color=2
;
; plot CRM model
;
;plot_io, d,ka(*,0), yra=[min(ka)>1,max(ka)<1e10],XTICKFORMAT='LABEL_DATE',xminor=12,$
;xtitle='UTC Date, Day of Year',xmargin=[10,4],ymargin=[-14.5,0],ytit='Proton Flux (p/s-cm2-sr-MeV)'
;if nreg1 gt 0 then oplot, de(reg1),ka(reg1),psym=3,color=4
;if nreg2 gt 0 then oplot, de(reg2),ka(reg2),psym=3,color=6
;if nreg3 gt 0 then oplot, de(reg3),ka(reg3),psym=3,color=5
;oplot, d,ka(*,0),color=2
;for i=0,n_elements(mot)-1 do oplot, [mot(i),mot(i)],10^!y.crange
;oplot, [d(0),d(n_elements(d)-1)],[5e4,5e4],color=1
;
;out=tvrd()
;WRITE_GIF, crmdir+'crmplxatt.gif', out,red,green,blue

end
