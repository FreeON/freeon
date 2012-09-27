---
layout: default
title: How To Flash the BIOS
---

We will boot from a bootable USB stick and then flash the BIOS with the appropriate flash utility.

BIOS Updates
------------

-   [AMD servers, TYAN FT48B8812 (B8812F48W8HR)](http://www.tyan.com/support_download_bios.aspx?model=B.FT48B8812): [Version 1.06 released on 2011/12/13](ftp://ftp.tyan.com/bios/FT48-B8812_v106.rar)

Create a bootable USB stick
---------------------------

First we need to create a bootable USB stick running MS-DOS.

-   [Rufus](http://rufus.akeo.ie/), which runs under Windows.
-   [UNetbootin](http://unetbootin.sourceforge.net/), which runs under windows or linux.
-   [Tutorial](http://www.sevenforums.com/tutorials/46707-ms-dos-bootable-flash-drive-create.html):
    -   Download the HP Flash Utility [hpflash1.zip](http://www.sevenforums.com/attachments/tutorials/42022d1260810265-ms-dos-bootable-flash-drive-create-hpflash1.zip) and also download the Windows 98 MS-DOS System Files [win98boot.zip](http://www.sevenforums.com/attachments/tutorials/42023d1260810265-ms-dos-bootable-flash-drive-create-win98boot.zip).
    -   Run the utility and select create MS-DOS disk pointing to the files extracted from the second zip file (see screenshot).

Install [FreeDOS](http://www.freedos.org/).
