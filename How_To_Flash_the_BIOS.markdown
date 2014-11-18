---
layout: default
title: How To Flash the BIOS
---

We will boot from a bootable USB stick and then flash the BIOS with the appropriate flash utility.

BIOS Updates
------------

-   Our 10 AMD nodes are [TYAN FT48B8812 (B8812F48W8HR)](http://www.tyan.com/product_SKU_spec.aspx?ProductType=BB&pid=434&SKU=600000186). The BIOS is available at [FT48B8812 BIOS](http://www.tyan.com/support_download_bios.aspx?model=B.FT48B8812): [Version 1.06 released on 2011/12/13](ftp://ftp.tyan.com/bios/FT48-B8812_v106.rar)

Creating a bootable USB stick
-----------------------------

First we need to create a bootable USB stick running MS-DOS.

![ thumb | Screenshot of the HP USB Disk Utility](HP USB Utility.png  " thumb | Screenshot of the HP USB Disk Utility")

-   [Rufus](http://rufus.akeo.ie/), which runs under Windows.
-   [UNetbootin](http://unetbootin.sourceforge.net/), which runs under Windows or Linux.
-   [Tutorial](http://www.sevenforums.com/tutorials/46707-ms-dos-bootable-flash-drive-create.html):
    -   Download the HP Flash Utility [hpflash1.zip](http://www.sevenforums.com/attachments/tutorials/42022d1260810265-ms-dos-bootable-flash-drive-create-hpflash1.zip) and also download the Windows 98 MS-DOS System Files [win98boot.zip](http://www.sevenforums.com/attachments/tutorials/42023d1260810265-ms-dos-bootable-flash-drive-create-win98boot.zip).
    -   Run the utility and select **"Create a DOS startup disk"** pointing to the files extracted from the second zip file (see screenshot).

-   [Another tutorial](http://www.chavers.us/robs-place-mainmenu-42/17-ubuntu-notes/46-easiest-way-to-create-a-usb-dos-boot-disk-using-linux):
    -   Partition USB stick.
    -   Install FreeBSD on the stick using UNetbootin under linux.
    -   Copy BIOS files into root of USB stick, will show up in drive B: or C: once running FreeDOS.

-   [Yet another tutorial](http://honeypot.net/2011/10/11/making-dos-usb-images-on-a-mac/):

So far USB either doesn't boot or I can't figure out how to add the BIOS files to the DOS live session.

Install [FreeDOS](http://www.freedos.org/).
