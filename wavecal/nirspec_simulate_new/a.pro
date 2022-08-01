;**********************************************************
; NAME:                      NIRSPEC_EFS
;
; PURPOSE:
; NIRSPEC_EFS.PRO
;     This program is the echelle format simulator for NIRSPEC.
;     It allows the user to setup a sequence of instrument
;     configurations and observing sequences.  Once setup, these
;     may either be saved as csh scripts, or executed immediately
;     at the telescope.  As a script generator, the EFS is ideal
;     for observers preparing for an observing run.  Scripts can
;     be sourced from any csh prompt on the instrument control
;     computer, or they can be reloaded into the EFS for
;     modifications. The csh scripts execute a series of "show"
;     and "modify" commands to communicate with the NIRSPEC server.
;     For directly driving the instrument, the EFS provides a "GO"
;     button that saves the current setup as a script in the standard
;     script directory and as part of number sequence of scripts.
;     The EFS then spawns a UNIX command that sources the script.
;     An "ABORT" button is also provided that initiates a different
;     csh command that kills scripting processes.
;
;     The EFS can run as a standalone package for script generation
;     or concurrently with the NIRSPEC server to drive the instrument.
;     Once a script is started, it loses all connection with the script.
;     In order to monitor the script, the user should use the
;     DATAVIEWS package and the Quicklook to view images.
;
;
; MODIFICATION HISTORY:
;
;       NIRSPEC_EFS started as a modified version of ECHWRAP_011 which
;	was a demo of Don Figer and Harry Teplitz's proposed redesign
;       to incorporate more User Interface features into to the Echelle
;       Format Simulator.  ECHWRAP interfaced to James Graham's
;	code for calculating echellograms.
;
;	version echgen_001: started Nov 12, 1997 by James Larkin/UCLA
;			    modified plotting algorithm to work in
;			    echelle and grating angles. Added many
;			    comments and changed gui placement and options.
;			    Still a stripped down version compared to
;			    the full operations of ECHWRAP.
;       nirspec_efs: Jan and Feb, 1998 by J. Larkin. Changed name as
;                           part of a new naming convention in which
;                           useful IDL scripts start with "nirspec"
;                           Scripting was changed to csh scripting.
;                           A starter program called run_efs was
;                           written to launch the idl in server mode
;                           and execute the EFS startup properly.
;                           Filters and slits were changed to pulldown
;                           menus.  Script reading and error checking
;                           were added, along with a number of
;                           minor modifications.
;       jim_efs:  Dec 2006 by J. Lyke.  added hooks for LGS-AO. mods
;                           start on line 1568 or so.
;       nod_efs:  Feb 2007 by J. Lyke.  Allows user-defined nod
;                          patterns
;       test_efs: Feb 2009 by J. Lyke.  Executes scripts in pop-up
;                          xterm
;       test_efs: Nov 2009 by J. Lyke. Change filter calls to use
;                          "filt" and "filter" commands to allow
;                          use of new AO pupil
;
;  
;       nirspec_efs: Dec 11, 2012 by G. Doppmann
;
;       Added some code
;       that had been removed after June 2011 to (re)-enable saving
;       setup scripts with the correct lines so that if could be read
;       back in by the EFS.
;
;
;       nirspec_efs: Dec 17, 2012 by G. Doppmann
;
;       Enabled setting keyword "script" to 1 during script execution
;        to safeguard against launching subsequent scripts.
;
;       nirspec_efs: Sept 17, 2014 by G. Doppmann
;
;       Disabled opening the hatch during setup mode.
;       This is to minmize dust getting onto the cal optics
;        before going on sky.  Now the hatch opens during
;        "Enable Nighttime Mode"
;
;
;       nirspec_efs:  Nov. 2018 GDoppmann
;                     Updated for H2RG detector params and new plate scale
;                     Display Works with HIRES mode, needs tweeks for lowres
;
;
;**********************************************************
; Important variables and their meanings:
;	nodmode:	specifies number of slit positions to
;			be used during exposures.  Nodmode
;                       is global to all setups, including the
;                       reference star frames.
;			values: 1 -> stare
;				2 -> 2 position nod
;				3 -> 3 position nod
;				4 -> 4 position nod
;                               5 -> ABBA
;				6 -> nod off to sky
;				7 -> user defined via a file
;
;       boxes is the array of structures that holds the configuration
;                       sequence information.
;
;	boxes(i).spmode:specifies spectroscopic mode.
;			values: 'high res' -> echelle w/ cross
;				'low res'  -> grating only
;				'imaging'  -> imaging mode (oth order)
;
;	calmode:	specifies if a full set of calibrations should
;			be done.
;			values: 0 -> Setup mode only
;                               1 -> Lamps only
;                               2 -> basic mode, only one object,
;				     no lamps.
;				3 -> Object and lamps only, no star
;                               4 -> quick mode, object and reference
;				    star, no lamps
;				5 -> full sequence: object, star, lamps
;
;
;	lambda:		array of 1000 locations along wavelength
;			range for the current filter
;
;	echphi:		array of 1000 angles for lambdas for the 
;			echelle = theta+delta
;
;	graphi:		array of 1000 angles for lambdas for the
;			grating = theta+deltas = alphas - ha
;

;	aspect:		the graphi/echphi aspect ratio for each filter
;			that makes boxes about square.
;
; Echellogram variables:
;	leftlam:	array of leftmost wavelengths for plotted orders
;	rightlam:	array of rightmost wavelengths for plotted orders.
;	leftech:	array of echelle angles for leftmost wavelengths.
;	rightech:	array of echelle angles for rightmost wavelengths.
;	leftgra:	array of grating angles for leftmost wavelengths.
;	rightgra:	array of grating angles for rightmost wavelengths.
;	numorders:	number of orders for current filter
;	uleftx:		array of echelle angles for upper slit ends at left.
;	ulefty:		array of grating angles for upper slit ends at left.
;	urightx:	array of echelle angles for upper slit ends at right.
;	urighty:	array of grating angles for upper slit ends at right.
;	lleftx:		array of echelle angles for lower slit ends at left.
;	llefty:		array of grating angles for lower slit ends at left.
;	lrightx:	array of echelle angles for lower slit ends at right.
;	lrighty:	array of grating angles for lower slit ends at right.
;
; Lowres variables:
;       botlam:         wavelength at bottom of screen
;       toplam:         wavelength at top of screen
;**********************************************************


;***********************************************************
; The following two functions are used for the rounding of 
; numbers to two decimal places.  This is particularly used
; for changes in the cross dispersor and the echelle grating, 
; which can only be moved to a precision of 0.01 degrees.  
; The first function rounds a number to the nearest hundredth, 
; the second converts an angle in radians to degrees, rounds
; to the nearest hundredth, then converts back to radians.
;***********************************************************

function round2, num

	num=round(num*100)/100.0
	return, num

end

function roundphi, num

        num=(!pi/180)*(round2(num*(180/!pi)))
        return, num

end

function neflw, x, y

        x=float(x)
        y=float(y)
        w=(464.951/2.0)+(463.740/2.0)+(0.00724829*x)+ $
          (1.30460e-05*x*x)+(0.0416329*y)+(5.00474e-06*y*y)
        return, w

end

function neflh, x, y

        x=float(x)
        y=float(y)
        h=(412.118/2.0)+(418.551/2.0)+(0.00173399*x)+ $
          (8.23022e-06*x*x)+(0.0168703*y)+(-1.44943e-05*y*y)
        return, h

end

function check_in_box, in_x, in_y, xs, ys, sptype

if sptype eq 'imaging' then begin
    x=in_x
    y=in_y
endif else begin
    temp=convert_coord(in_x, in_y, /device, /to_data)
    x=temp(0)
    y=temp(1)
endelse

xmax=(xs(2)+((xs(1)-xs(2))/(ys(1)-ys(2)))*(y-ys(2)))
xmin=(xs(3)+((xs(0)-xs(3))/(ys(0)-ys(3)))*(y-ys(3)))
ymax=(ys(3)+((ys(2)-ys(3))/(xs(2)-xs(3)))*(x-xs(3)))
ymin=(ys(0)+((ys(1)-ys(0))/(xs(1)-xs(0)))*(x-xs(0)))
if (x lt xmax and x gt xmin and y lt ymax and y gt ymin) then $
  return, 1 $
else return, 0

end
	
;**********************************************************

;**********************************************************
PRO formatbase_event, Event
; A default routine that accepts and ignores all events that are
; outside of buttons and the draw box.   Since it is the first
; procedure the common blocks are defined within it.
;**********************************************************
common shared_1, font, echtheta, echsigma, dlambda, fullname, $
  nodmode, calmode, echgamma, gratheta, leftlam, rightlam, $
  leftech, rightech, leftgra, rightgra, numorders, numfil, filename, $
  userfile
common shared_2, totbox, boxes, data_x1, data_x2, order1, order2, $
  data_y1, data_y2, few, feh, grahaang,graorder, uleftx, ulefty, $
  urightx, urighty, yarc, xarc, lleftx, llefty, lrightx, lrighty
common shared_3, boxnum, formatbase, formatmenu, slitmenu, thinmenu, $
  obj_name, star_name, itime, stime, coadds, filtermenu, yaspect, $
  pix, npix, maxbox, curbox, echphi, graphi, filterarr, botlam, toplam, $
  filtercmd
common shared_4, new_box, draw, full, nods, specmode, slitlen, $
  l1, l2, grasigma, moving, pre_box, next_box, aspect, del_box, $
  sel_col, off_col, spec_col, border_col, text_col, status_text, $
  echbox, grabox, lambox, ordbox, xpixbox, ypixbox, sfilenum, scriptname, $
  scriptdir, scr_name_base, scriptdir2
common shared_5, slitlist, save_button, scoadds, chipang, num_scr, $
  plot_x1, plot_x2, plot_y1, plot_y2, lgraphi, hgraphi, scrx, scry, serv, $
  dist_few, dist_feh, reps, oh_lines, oh_stat, neon_lines, neon_stat, $
  krypton_lines, krypton_stat, argon_lines, argon_stat, xenon_lines, $
  xenon_stat, private_lines, private_stat, oh_num, neon_num, krypton_num, $
  argon_num, xenon_num, private_num, oh_col, argon_col, xenon_col, $
  krypton_col, neon_col, user_col, plot_filter_button, filtermenu_index, $
  user_stat, user_list_length, box_button_base, bad, numbad, mask_flag, $
  mask_toggle, i_menu, blank, xoff, yoff, box_x, box_y, thinname, i_boxbase, $
  valid_boxes, current_fil, sptype, current_slit, curbox_field, csr_x, csr_y, $
  in_draw, working_dir, lock_toggle, lock_flag
common shared_6, filelist, new_line, filelist_box, add_line_box, $
		line_list, data_list, save_list_box, current_dir, save_file, $ 
		save_list_name, color_base, rgb, r, g, b, color_id, draw_id, $
		browse_button, old_list, result, displayed_list, comments, $
		add_comment_box, old_comments, list, user_lines_id, z_box, $
                shifted_list, redshift, changed, user_base
common plot_shared, plotted, efs_base, efs_draw_index
common shared_go, go_button
common shared_ao, aomode, targwave, string_tw

; if a top level base event, then do this routine
if event.id eq event.top then begin
    ; Incase this is a resize event, get the event coordinates. JLW
    base_x=event.x
    base_y=event.y

    ; Check if the size is less than 608x686 pixels
    if base_x lt 608 then base_x=608
    if base_y lt 686 then base_y=686

    ; The graphics window is smaller than the overall window
    scrx=base_x-256
    scry=base_y-312

    ; Set the draw window's size
    widget_control, draw, xsize=scrx
    widget_control, draw, ysize=scry

    ; Set the overall size
    widget_control, event.top, xsize=base_x-4
    widget_control, event.top, ysize=base_y-40

    ; Redraw the echellogram
    draw_echellogram

    ; Draw the current set of array positions
    draw_boxes

    ; call fill_curbox_arrays to update box_x and box_y with the new
    ; data/device coordinate conversion
    fill_curbox_arrays

endif

END
;**********************************************************

;**********************************************************
PRO formatmenu_event, Event
; This procedure accepts all of the main menu events, and sends
; them to the appropriate routines.
;**********************************************************
common shared_1
common shared_2
common shared_3
common shared_4
common shared_5
common shared_6


widget_control, Event.id, get_value=val
CASE val OF 
; File menu items
    'Open Configuration...     <Ctrl+o>': open_config
    'Save Configuration As... <Ctrl+s>': save_button_event
    'Configure Script File       <Ctrl+f>': config_script
    'Quit                                       <Ctrl+q>': format_quit
 
; Telescope menu items JLW
    'Set Nod Parameters  <Ctrl+p>':begin
        ;print, event.id
        print, 'run_nod'
        spawn, 'run_nod &'
        end

; Overlays
    'Clear Overlays  <Ctrl+c>': begin
        ;print, event.id
        
        ; Clear checks on arc line menus  
        widget_control, event.id+4, set_value='[  ]   Neon         <Ctrl+n>'
        widget_control, event.id+5, set_value='[  ]   Argon       <Ctrl+a>'
        widget_control, event.id+6, set_value='[  ]   Krypton  <Ctrl+k>'
        widget_control, event.id+7, set_value='[  ]   Xenon       <Ctrl+x>'
        widget_control, event.id+2, set_value='[  ]   OH Lines    <Ctrl+l>'
        widget_control, event.id+1, set_value='[  ]   User Lines  <Ctrl+u>'
        widget_control, i_menu+6 , set_value='user'
        widget_control, i_menu+7 , set_value='oh'
        widget_control, i_menu+8 , set_value='neon'
        widget_control, i_menu+9 , set_value='argon'
        widget_control, i_menu+10 , set_value='krypton'
        widget_control, i_menu+11 , set_value='xenon'
        clear_overlay
        end
    '[  ]   User Lines  <Ctrl+u>': begin
        widget_control, event.id, set_value='[x]   User Lines  <Ctrl+u>'
        widget_control, i_menu+6 , set_value='user2'
        user_lines_id=event.id ;save id number for cancel button
        user_lines_setup
        end
    '[x]   User Lines  <Ctrl+u>': begin
        widget_control, event.id, set_value='[  ]   User Lines  <Ctrl+u>'
        widget_control, i_menu+6 , set_value='user'
        user_lines_setup
        end
;    'Wavelength...': input_wavelength
    '[  ]   OH Lines    <Ctrl+l>': begin
        widget_control, event.id, set_value='[x]   OH Lines    <Ctrl+l>'
        widget_control, i_menu+7 , set_value='oh2'
        oh_lines_setup
        end
    '[x]   OH Lines    <Ctrl+l>': begin
        widget_control, event.id, set_value='[  ]   OH Lines    <Ctrl+l>'
        widget_control, i_menu+7 , set_value='oh'
        oh_lines_setup
        end
    '[  ]   Neon         <Ctrl+n>': begin
        ;print, event.id
        widget_control, event.id, set_value='[x]   Neon         <Ctrl+n>'
        widget_control, i_menu+8 , set_value='neon2'
        neon_lines_setup
        end
    '[x]   Neon         <Ctrl+n>': begin
        widget_control, event.id, set_value='[  ]   Neon         <Ctrl+n>'
        widget_control, i_menu+8 , set_value='neon'
        neon_lines_setup
        end
    '[  ]   Argon       <Ctrl+a>': begin
        widget_control, event.id, set_value='[x]   Argon       <Ctrl+a>'
        widget_control, i_menu+9 , set_value='argon2'
        argon_lines_setup
        end
    '[x]   Argon       <Ctrl+a>': begin
        widget_control, event.id, set_value='[  ]   Argon       <Ctrl+a>'
        widget_control, i_menu+9 , set_value='argon'
        argon_lines_setup
        end
    '[  ]   Krypton  <Ctrl+k>': begin
        widget_control, event.id, set_value='[x]   Krypton  <Ctrl+k>'
        widget_control, i_menu+10 , set_value='krypton2'
        krypton_lines_setup
        end
    '[x]   Krypton  <Ctrl+k>': begin
        widget_control, event.id, set_value='[  ]   Krypton  <Ctrl+k>'
        widget_control, i_menu+10 , set_value='krypton'
        krypton_lines_setup
        end
    '[  ]   Xenon       <Ctrl+x>': begin
        widget_control, event.id, set_value='[x]   Xenon       <Ctrl+x>'
        widget_control, i_menu+11 , set_value='xenon2'
        xenon_lines_setup
        end
    '[x]   Xenon       <Ctrl+x>': begin
        widget_control, event.id, set_value='[  ]   Xenon       <Ctrl+x>'
        widget_control, i_menu+11 , set_value='xenon'
        xenon_lines_setup
        end
; Help
    'Keyboard Hotkey List': begin
        hotkey_list=['NIRSPEC ECHELLE FORMAT SIMULATOR', $
                     '      KEYBOARD HOTKEY LIST', $
                     '', $
                     ' FILE', $
                     '     Open Configuration     Ctrl+o', $
                     '     Save Configuration     Ctrl+s', $
		     '     Configure Script File  Ctrl+f', $
                     '     Quit                   Ctrl+q', $
                     '', $
                     ' TELESCOPE', $
                     '     Set Nod Parameters     Ctrl+p', $
                     '', $
                     ' OVERLAYS', $
                     '     Clear Overlays         Ctrl+c', $
                     '     User Lines             Ctrl+u', $
                     '     OH Lines               Ctrl+l', $
                     '   ARC LINES', $
                     '         Neon               Ctrl+n', $
                     '         Argon              Ctrl+a', $
                     '         Krypton            Ctrl+k', $
                     '         Xenon              Ctrl+x', $
                     '', $
                     ' BUTTONS', $
                     '     Go                     Ctrl+g', $
                     '     Abort                  Meta+x', $
                     '   SETUP BOXES', $
                     '         New Box            Insert', $
                     '         Previous Box       Pageup', $
                     '         Next Box         PageDown', $
                     '         Delete Box         Delete', $
                     '', $
                     ' ARROW KEYS:  Use Arrow Keys to', $
                     '    move box by one pixel.', $
                     '']
        answer=dialog_message(hotkey_list, dialog_parent=formatmenu, $
                             title='EFS Hotkeys', /information)
        end
    'User Manual': BEGIN
        command='netscape help/efs_toc.html &'
        spawn, command
        end

ENDCASE

END
;**********************************************************
;**********************************************************
PRO keyboard_event, Event
; This procedure accepts all of the keyboard events, and sends
; them to the appropriate routines.
;**********************************************************
common shared_1
common shared_2
common shared_3
common shared_4
common shared_5
common shared_6

; if statement for arrow keys for what to pass to calc_pixel
if in_draw then begin
    pix_x=csr_x
    pix_y=csr_y
endif else begin
    pix_x=box_x(curbox)
    pix_y=box_y(curbox)
endelse

widget_control, Event.id, get_value=val
CASE val OF 
; File menu items
    'open': open_config
    'save': save_button_event
    'cfg_script' :config_script
    'quit': format_quit
 
; Telescope menu items JLW
    'nod':begin
        print, 'run_nod'
        spawn, 'run_nod &'
        end

; Overlays
    'clear': begin
        ;print, event.id
        ; Clear checks on arc line menus  
        widget_control, formatmenu+12, set_value='[  ]   Neon         <Ctrl+n>'
        widget_control, formatmenu+13, set_value='[  ]   Argon       <Ctrl+a>'
        widget_control, formatmenu+14, set_value='[  ]   Krypton  <Ctrl+k>'
        widget_control, formatmenu+15, set_value='[  ]   Xenon       <Ctrl+x>'
        widget_control, formatmenu+10, set_value='[  ]   OH Lines    <Ctrl+l>'
        widget_control, formatmenu+9, set_value='[  ]   User Lines  <Ctrl+u>'
        widget_control, event.id+3, set_value='neon'
        widget_control, event.id+4, set_value='argon'
        widget_control, event.id+5, set_value='krypton'
        widget_control, event.id+6, set_value='xenon'
        widget_control, event.id+2, set_value='oh'
        widget_control, event.id+1, set_value='user'
        clear_overlay
        end
    'user': begin
        widget_control, formatmenu+9, set_value='[x]   User Lines  <Ctrl+u>'
        widget_control, event.id, set_value='user2'
        user_lines_id=formatmenu+9 ;save id number for cancel button
        user_lines_setup
        end
    'user2': begin
        widget_control, formatmenu+9, set_value='[  ]   User Lines  <Ctrl+u>'
        widget_control, event.id, set_value='user'
        user_lines_setup
        end
;    'Wavelength...': input_wavelength
    'oh': begin
        widget_control, formatmenu+10, set_value='[x]   OH Lines    <Ctrl+l>'
        widget_control, event.id, set_value='oh2'
        oh_lines_setup
        end
    'oh2': begin
        widget_control, formatmenu+10, set_value='[  ]   OH Lines    <Ctrl+l>'
        widget_control, event.id, set_value='oh'
        oh_lines_setup
        end
    'neon': begin
        ;print, formatmenu
        widget_control, formatmenu+12, set_value='[x]   Neon         <Ctrl+n>'
        widget_control, event.id, set_value='neon2'
        neon_lines_setup
        end
    'neon2': begin
        widget_control, formatmenu+12, set_value='[  ]   Neon         <Ctrl+n>'
        widget_control, event.id, set_value='neon'
        neon_lines_setup
        end
    'argon': begin
        widget_control, formatmenu+13, set_value='[x]   Argon       <Ctrl+a>'
        widget_control, event.id, set_value='argon2'
        argon_lines_setup
        end
    'argon2': begin
        widget_control, formatmenu+13, set_value='[  ]   Argon       <Ctrl+a>'
        widget_control, event.id, set_value='argon'
        argon_lines_setup
        end
    'krypton': begin
        widget_control, formatmenu+14, set_value='[x]   Krypton  <Ctrl+k>'
        widget_control, event.id, set_value='krypton2'
        krypton_lines_setup
        end
    'krypton2': begin
        widget_control, formatmenu+14, set_value='[  ]   Krypton  <Ctrl+k>'
        widget_control, event.id, set_value='krypton'
        krypton_lines_setup
        end
    'xenon': begin
        widget_control, formatmenu+15, set_value='[x]   Xenon       <Ctrl+x>'
        widget_control, event.id, set_value='xenon2'
        xenon_lines_setup
        end
    'xenon2': begin
        widget_control, formatmenu+15, set_value='[  ]   Xenon       <Ctrl+x>'
        widget_control, event.id, set_value='xenon'
        xenon_lines_setup
        end
    'go': go_button_event, event
    'abort': abort_button_event, event
    'newbox': new_box_event, event
    'delbox': del_box_event, event
    'prebox': pre_box_event, event
    'nextbox': next_box_event, event

; Arrow Keys
    'up': begin
        calc_pixel, pix_x, pix_y
        boxes(curbox).graphi = boxes(curbox).graphi+(0.01*!pi/180)
        if in_draw then box_y(curbox)=box_y(curbox)+1
        erase_curbox            ; Erase current box
        find_corners
        erase_curbox            ; Draw current box
        if mask_flag then begin
            erase_bad
            find_bad
            erase_bad
        endif
        widget_control, grabox, set_value=(boxes(curbox).graphi+grahaang)*180/!pi
        end
    'down': begin
        calc_pixel, pix_x, pix_y
        boxes(curbox).graphi = boxes(curbox).graphi-(0.01*!pi/180)
        if in_draw then box_y(curbox)=box_y(curbox)-1
        erase_curbox            ; Erase current box
        find_corners
        erase_curbox            ; Draw current box
        if mask_flag then begin
            erase_bad
            find_bad
            erase_bad
        endif
        widget_control, grabox, set_value=(boxes(curbox).graphi+grahaang)*180/!pi
        end
    'left': begin
        calc_pixel, pix_x, pix_y
;        boxes(curbox).echphi = boxes(curbox).echphi-(0.01*!pi/180)
;        if in_draw then box_x(curbox)=box_x(curbox)-1
        if in_draw then begin
            boxes(curbox).echphi = boxes(curbox).echphi-(0.01*!pi/180)
            box_x(curbox)=box_x(curbox)-1
            erase_curbox        ; Erase current box
            find_corners
            erase_curbox        ; Draw current box
            if mask_flag then begin
                erase_bad
                find_bad
                erase_bad
            endif
            widget_control, echbox, set_value=boxes(curbox).echphi*180/!pi
        endif
    end
    'right': begin
        calc_pixel, pix_x, pix_y
;        boxes(curbox).echphi = boxes(curbox).echphi+(0.01*!pi/180)
;        if in_draw then box_x(curbox)=box_x(curbox)+1
        if in_draw then begin
            boxes(curbox).echphi = boxes(curbox).echphi+(0.01*!pi/180)
            box_x(curbox)=box_x(curbox)+1
            erase_curbox        ; Erase current box
            find_corners
            erase_curbox        ; Draw current box
            if mask_flag then begin
                erase_bad
                find_bad
                erase_bad
            endif
        endif
        widget_control, echbox, set_value=boxes(curbox).echphi*180/!pi
        end

ENDCASE

END
;**********************************************************

;**********************************************************
PRO filtermenu_event, Event, findex=findex
; This procedure handles the pulldown menu events for the
; filter wheels.
; It just takes the number of the menu pulled, and takes
; that filter entry and places it in the boxes.filter
; structure location.  Then it updates the standard procedures
; for updating the status of the buttons, selecting filter
; for the echellogram, calculate a new echellogram, draw
; the active array on the echellogram and any others that
; share that filter.
;**********************************************************

common shared_1
common shared_2
common shared_3
common shared_4
common shared_5
common shared_6
common shared_ao

widget_control, draw, get_value=index
wset, index

filtermenu_index(curbox)=Event.index
current_fil=event.index

help, current_fil

; Set blocking mode to AO stop, or thin
; jlyke 2009 nov 17 don't think we have to change this...
;if (aomode) then begin
; boxes(curbox).blocker=2 
;endif else begin
; widget_control, thinmenu, set_value=1
; boxes(curbox).blocker=1
;endelse

; Set the boxes filter and wheel parameters
;boxes(curbox).filter = filterarr(Event.index).position
;boxes(curbox).wheel = filterarr(Event.index).wheel
boxes(curbox).filter = filterarr(current_fil).position
boxes(curbox).wheel = filterarr(current_fil).wheel

; Determine targwave for DAR correction
; do not send it until user presses "GO"

;l1 = filterarr(Event.index).l1
;l2 = filterarr(Event.index).l2
l1 = filterarr(current_fil).l1
l2 = filterarr(current_fil).l2
; wavelengths in filterarr are in mm, targwave is in micron
targwave = 1000.*(l1 + l2)/2.
string_tw = strcompress(string(targwave), /remove_all)
widget_control, set_value = 'Targwave will be '+string_tw, status_text


blk=where(filterarr.name eq 'BLANK')
;if event.index eq blk(0) then begin  ; blank selected
if current_fil eq blk(0) then begin  ; blank selected
    blank=1

    ; disable buttons
    widget_control, slitmenu, sensitive=0
    widget_control, echbox, sensitive=0
    widget_control, grabox, sensitive=0
    widget_control, plot_filter_button, sensitive=0
	
    widget_control, lambox, set_value='Undispersed'
    widget_control, ordbox, set_value='None'
    widget_control, xpixbox, set_value='N/A'
    widget_control, ypixbox, set_value='N/A'


    if mask_flag ne 0 then begin
        mask_flag=0
        widget_control, mask_toggle, set_value=1
        erase_bad
    endif

endif else begin     ; not blank
    blank=0

    widget_control, slitmenu, sensitive=1
    widget_control, echbox, sensitive=1
    widget_control, grabox, sensitive=1

    ; Update the buttons and calculate the echellogram.
    update_buttons
    choose_standard_filter
    setup_ech_filter

; Center current box in echellogram
if ( boxes(curbox).spmode eq 'high res' ) then begin
    boxes(curbox).echphi = roundphi(echtheta)
    boxes(curbox).graphi = roundphi(graphi(499))
endif
if ( boxes(curbox).spmode eq 'low res' ) then begin
    boxes(curbox).graphi = roundphi(graphi(499))
; Changed temporarily for first light to 179.99 degrees
;    boxes(curbox).echphi = roundphi(180*!dtor)
    boxes(curbox).echphi = roundphi(179.7*!dtor)
endif
if ( boxes(curbox).spmode eq 'imaging' ) then begin
    boxes(curbox).graphi = roundphi(0.0*!pi/180.0)
; Changed temporarily for first light to 179.99 degrees
;    boxes(curbox).echphi = roundphi(!pi)
    boxes(curbox).echphi = roundphi(179.7*!dtor)
endif

find_corners

endelse

; Draw the echellogram and all setups in this filter.

fill_curbox_arrays
calc_pixel, box_x(curbox), box_y(curbox)
draw_echellogram
draw_boxes


END
;**********************************************************

;**********************************************************
PRO fill_curbox_arrays
; added by jlw 2/11/99.  called when a new config is opened
;**********************************************************
common shared_1
common shared_2
common shared_3
common shared_4
common shared_5
common shared_6
common shared_ao

for n=1, totbox do begin
        ; get filter
	for i = 0, (numfil-1) do $
		if ( (boxes(n).filter eq filterarr(i).position) and $
		(boxes(n).wheel eq filterarr(i).wheel) ) then $
			filtermenu_index(n)=i
        if boxes(n).spmode ne 'imaging' and filtermenu_index(n) ne 18 $
          then begin
            ; enter center of box in data coords
            tmp=[(boxes(n).x(0)+boxes(n).x(2))/2.0, $
                 (boxes(n).y(0)+boxes(n).y(2))/2.0]
            widget_control, draw, get_value=index
            wset, index
            temp=convert_coord(tmp(0), tmp(1), /data, /to_device)
            box_x(n)=temp(0)
            box_y(n)=temp(1)
        endif else begin
            if box_x(n) le 0 or box_y(n) le 0 or box_x(n) ge scrx or box_y(n) $
              gt (scry) then begin
                box_x(n)=scrx/2.0
                box_y(n)=scry/2.0
                draw_echellogram
                draw_boxes
            endif
        endelse
endfor

end
;**********************************************************

;**********************************************************
PRO thinmenu_event, Event
; Handles changes to the blocking menu.  Updates the value
; to boxes(curbox).blocker.
;**********************************************************
common shared_1
common shared_2
common shared_3
common shared_4
common shared_5
common shared_6
common shared_ao

; get new selection and fill it into boxes structure

if (aomode) then begin
 boxes(curbox).blocker=blocker 
endif else begin
 widget_control, thinmenu, get_value=blocker
 boxes(curbox).blocker=blocker
; filtermenu_event, event, findex=current_fil
endelse

end

;**********************************************************
 PRO plot_filter_button_event, Event
; Calls exeternal plot_filter procedure which displays the 
; profile of the selected filter's throughput.
;**********************************************************
common shared_1
common shared_2
common shared_3
common shared_4
common shared_5
common shared_6
common shared_ao
 
; gather info about echelle and cross dispersor
echinfo={sigma:echsigma, gamma:echgamma, theta:echtheta}
grainfo={sigma:grasigma, theta:gratheta, haang:grahaang}
; pass to plot_filter to do the plotting
;temp_plot_filter, current_fil, echinfo, grainfo
;plot_filter, current_fil, echinfo, grainfo
plot_filter, filterarr(current_fil), echinfo, grainfo

end

;**********************************************************
 PRO lock_toggle_event, Event
; Determines what to do if lock toggle is hit.
; event structure is {ID:0L, TOP:0L, HANDLER:0L, VALUE:0}
; where VALUE is either 0 (unlock) or 1 (lock), 
; the index of the button pressed.
;**********************************************************
common shared_1
common shared_2
common shared_3
common shared_4
common shared_5
common shared_6
common shared_ao

; print, event.value

if (event.value eq 0 ) then begin
  ; unlock
  lock_flag = 0
  widget_control, specmode, sensitive=1
  widget_control, slitmenu, sensitive=1
  widget_control, filtermenu, sensitive=1
  if (not aomode) then begin
    widget_control, thinmenu, sensitive=1
  endif
  widget_control, echbox, sensitive=1
  widget_control, grabox, sensitive=1
  widget_control, draw, sensitive=1
endif else begin
  ; lock
  lock_flag = 1
  widget_control, specmode, sensitive=0
  widget_control, slitmenu, sensitive=0
  widget_control, filtermenu, sensitive=0
  if (not aomode) then begin
    widget_control, thinmenu, sensitive=0
  endif
  widget_control, echbox, sensitive=0
  widget_control, grabox, sensitive=0
  widget_control, draw, sensitive=0
endelse

end

;**********************************************************
 PRO mask_toggle_event, Event
; Determines what to do if mask toggle is hit.
; event structure is {ID:0L, TOP:0L, HANDLER:0L, VALUE:0}
; where VALUE is either 0 or 1, the index of the button 
; pressed.
;**********************************************************
common shared_1
common shared_2
common shared_3
common shared_4
common shared_5
common shared_6

; print, event.value

if event.value eq 0 then begin
    mask_flag=1
    find_bad
    erase_bad
endif else begin
    mask_flag=0
    erase_bad
endelse

end

pro find_bad
     common shared_1
     common shared_2
     common shared_3
     common shared_4
     common shared_5
     common shared_6


; print, boxes(curbox).spmode
bechphi=boxes(curbox).echphi
bgraphi=boxes(curbox).graphi
for i=0, numbad-1 do begin
    if boxes(curbox).spmode eq 'high res' then begin
        bad(i).plot_x=asin((sin(bechphi)+sin(bechphi+bad(i).echdel))/2.0)
        bad(i).plot_y=asin((sin(bgraphi+grahaang)+sin(bgraphi $
             -grahaang+bad(i).gradel))/(2.0*cos(grahaang)))
    endif 
    if boxes(curbox).spmode eq 'low res' then begin
        bad(i).plot_x=bad(i).x+200.0
        bad(i).plot_y=asin((sin(bgraphi+grahaang)+sin(bgraphi $
             -grahaang+bad(i).gradel))/(2.0*cos(grahaang)))
    endif
endfor    

end

pro erase_bad
     common shared_1
     common shared_2
     common shared_3
     common shared_4
     common shared_5
     common shared_6

widget_control, draw, get_value=index
wset, index

if ( boxes(curbox).spmode ne 'imaging' ) then begin
    device, set_graphics = 6
    col = sel_col
    plots, bad.plot_x, bad.plot_y, color=col, psym=3
    device, set_graphics=3
endif

end


;**********************************************************
PRO nods_event, Event
; Handles changes to the nod selector.  Nod controls the
; number of exposures/positions along the slit to use for
; each grating setup.  It is global in that all setups will
; share this nodmode.
;**********************************************************
common shared_1
common shared_2
common shared_3
common shared_4
common shared_5
common shared_6
common shared_ao

; get new value and save to nodmode variable
print, '****************************'
nodmode = event.index + 1
print, 'nod mode =', nodmode
widget_control, set_value = "Ready to set up instrument", status_text

if ( nodmode eq 7 ) then begin
  userfile = dialog_pickfile( $
    path=working_dir, get_path=new_dir, $
    title='Select user-defined nod pattern file', $
    /must_exist )
  tmp = modify('nirspec.userstr3', userfile)
  j=strlen(userfile)
  if (j gt 0) then begin
    working_dir=new_dir
    widget_control, set_value = "Nod file: "+userfile, status_text
  endif else begin
    widget_control, set_value = "No nod file selected", status_text
  endelse
endif

END
;**********************************************************


;**********************************************************
PRO reps_event, Event
; Handles changes to the repitions selector.  Repetitions controls the
; times a nod pattern is repeated at each setting.  It can be set for
; each box setup separately.
;**********************************************************
common shared_1
common shared_2
common shared_3
common shared_4
common shared_5
common shared_6
common shared_ao

; get new value and save to boxes.reps variable
print, '****************************'
boxes(curbox).reps = event.index + 1
boxes(0).reps=event.index+1
print, 'rep mode =', boxes(curbox).reps

END
;**********************************************************

;**********************************************************
PRO full_event, Event
; Handles changes to the calibration mode selector.  The
; calmode variable determines if a reference star is observed
; as part of the script and if automatic calibration frames
; are taken after the observations.  The available running 
; options, such as setting parameters like star names, coadds, 
; etc. depend on which mode is selected, and therefore, some
; options are disabled for certain modes.
;**********************************************************
common shared_1
common shared_2
common shared_3
common shared_4
common shared_5
common shared_6
common shared_ao


print, '****************************'
calmode = event.index

case calmode of 
    0:begin                     ; Setup Only
        widget_control, nods, sensitive=0
        widget_control, reps, sensitive=0
        widget_control, obj_name, sensitive=0
        widget_control, star_name, sensitive=0
        widget_control, itime, sensitive=0
        widget_control, stime, sensitive=0
        widget_control, coadds, sensitive=0
        widget_control, scoadds, sensitive=0
;        widget_control, boxnum, sensitive=0
        widget_control, box_button_base, sensitive=0
        widget_control, i_boxbase, sensitive=0
        curbox=1
        totbox=1
        widget_control, curbox_field, set_value=strcompress(1)
        update_buttons
        calc_pixel, box_x(curbox), box_y(curbox)
        choose_standard_filter
        setup_ech_filter
        if mask_flag then erase_bad
        find_corners
        draw_echellogram
        draw_boxes
    end
    1:begin                     ; Lamps Only
        widget_control, nods, sensitive=0
        widget_control, reps, sensitive=0
        widget_control, obj_name, sensitive=0
        widget_control, star_name, sensitive=0
        widget_control, itime, sensitive=0
        widget_control, stime, sensitive=0
        widget_control, coadds, sensitive=0
        widget_control, scoadds, sensitive=0
;        widget_control, boxnum, sensitive=0
        widget_control, box_button_base, sensitive=0
        widget_control, i_boxbase, sensitive=0
        curbox=1
        totbox=1
        widget_control, curbox_field, set_value=strcompress(1)
        update_buttons
        calc_pixel, box_x(curbox), box_y(curbox)
        choose_standard_filter
        setup_ech_filter
        if mask_flag then erase_bad
        find_corners
        draw_echellogram
        draw_boxes
    end
    2:begin                     ; Object Only
        widget_control, nods, sensitive=1
        widget_control, reps, sensitive=1
        widget_control, obj_name, sensitive=1
        widget_control, star_name, sensitive=0
        widget_control, itime, sensitive=1
        widget_control, stime, sensitive=0
        widget_control, coadds, sensitive=1
        widget_control, scoadds, sensitive=0
;        widget_control, boxnum, sensitive=1
        widget_control, box_button_base, sensitive=1
        widget_control, i_boxbase, sensitive=1
        update_buttons
    end
    3:begin                     ; Bright Object Only
        widget_control, nods, sensitive=1
        widget_control, reps, sensitive=1
        widget_control, obj_name, sensitive=1
        widget_control, star_name, sensitive=0
        widget_control, itime, sensitive=1
        widget_control, stime, sensitive=0
        widget_control, coadds, sensitive=1
        widget_control, scoadds, sensitive=0
;        widget_control, boxnum, sensitive=1
        widget_control, box_button_base, sensitive=1
        widget_control, i_boxbase, sensitive=1
        update_buttons
    end
endcase

print, 'cal mode =', calmode

END

;**********************************************************

;**********************************************************
PRO slitmenu_event, Event
common shared_1
common shared_2
common shared_3
common shared_4
common shared_5
common shared_6
common shared_ao

; get new slits
current_slit=event.index
boxes(curbox).slit = Event.index
boxes(0).slit=event.index
slitlen = slitlist(boxes(curbox).slit).length

; update gui
update_buttons
setup_ech_filter

; find where new box should be drawn and draw it
find_corners
draw_echellogram
draw_boxes


END	
;**********************************************************

;**********************************************************
PRO specmode_event, Event
common shared_1
common shared_2
common shared_3
common shared_4
common shared_5
common shared_6
common shared_ao


; get spectral mode and update gui and echellogram accordingly
spmode = Event.index

case spmode of
	0:sptype='high res'
	1:sptype='low res'
	2:sptype='imaging'
endcase

boxes(curbox).spmode=sptype

if spmode ne 2 then begin
; Changed temporarily for first light to 179.99 degrees
    boxes(curbox).echphi = roundphi(179.7*!dtor)
;    boxes(curbox).echphi=roundphi(180*!dtor)
    if spmode eq 0 then boxes(curbox).echphi=roundphi(echtheta)
    boxes(curbox).graphi=roundphi(graphi(499))
    widget_control, mask_toggle, sensitive=1
    calc_pixel, box_x(curbox), box_y(curbox)
    update_buttons
    choose_standard_filter
    setup_ech_filter
    find_corners
    draw_echellogram
    draw_boxes
endif else begin
; Changed temporarily for first light to 179.99 degrees
    boxes(curbox).echphi = roundphi(179.7*!dtor)
;    boxes(curbox).echphi = roundphi(180*!dtor)
    boxes(curbox).graphi = roundphi(0.0)
;    widget_control, draw, get_value=index
;    wset, index
;     widget_control, mask_toggle, sensitive=0
;    plot, [1,46], [1,46], xtitle='arcseconds', $
;      ytitle='arcseconds', xstyle=1, ystyle=1, /nodata
;    xyouts, 23, 23, 'IMAGING', /data, charsize=8, align=0.5, $
;      charthick=10
    calc_pixel, box_x(curbox), box_y(curbox)
    update_buttons
    draw_echellogram
    draw_boxes
endelse

END

;**********************************************************


;**********************************************************
PRO abort_button_event, Event
; Abort button has been pushed.  Activate the abort_script
; csh routine that kills script processes.

common shared_1
common shared_2
common shared_3
common shared_4
common shared_5
common shared_6
common shared_ao

inp = strarr(2)


print, 'aborting script'

; Now spawn to execute the command
inp(0) = 'abort_script &'
msg = 'Aborting script.'
widget_control, set_value = msg, status_text
spawn, inp(0)

end

;**********************************************************
PRO go_button_event, Event
common shared_1
common shared_2
common shared_3
common shared_4
common shared_5
common shared_6
common shared_ao

print, 'Selecting GO'
widget_control, set_value = "You've selected to go.", status_text
num_scr = num_scr + 1
date = string(systime())        ; Get the system date.
dd = strmid(date, 8, 2)
mm = strmid(date, 4, 3)
scriptname = 's' +  string(format='(I2.2)',dd) + $
  mm + string(format='(I4.4)',num_scr) + '.csh'
;inp(1) = '/kroot/kss/nirspec/ui/xnirspec/scripts/' + scriptname
;scriptdir = '/kroot/kss/nirspec/ui/xnirspec/scripts'
inp = strarr(2)
;inp(1) = scriptdir + '/' + scriptname

; MODIFIED   now writes to the current directory  GWD
;inp(1) = scriptdir + '/' + scriptname
inp(1) = scriptname

on_ioerror, bad
unit = 0
; Write out the script file in order:
filename = inp(1)
;sfilenum = show("nirspec.filenum", /ascii)
;i = ktlClose("nirspec")
openw, unit, filename, /get_lun
save_script, unit
;close, unit
free_lun, unit

; Now spawn to execute the command
inp(0) = 'xterm -e sleep_wrapper -n 1 source ' + inp(1) + ' &'
print, inp(0)
spawn, inp(0)
inp(0) = 'Starting script ' + inp(1) + '.'
widget_control, set_value = inp(0), status_text

return
bad: on_ioerror, null
if ( unit gt 0 ) then free_lun, unit
err = dialog_message("IO Error has occurred writing script.", $
                     dialog_parent=formatbase)
print, "IO Error has occured writing in script."
widget_control, set_value = "Error writing out script.", status_text

END
;**********************************************************

;**********************************************************
PRO save_button_event, Event
common shared_1
common shared_2
common shared_3
common shared_4
common shared_5
common shared_6
common shared_ao

on_ioerror, bad
unit = 0
new_path=''
filename = dialog_pickfile(filter='*', get_path=new_path, path=working_dir)
fdecomp, filename, sv_disk, sv_dir, sv_name, sv_ext  
;scriptdir = '.'
scriptname = strcompress(sv_name+'.'+sv_ext, /remove_all)
j = strlen(filename)
if ( j gt 0 ) then begin
    working_dir=new_path
    scriptdir2=new_path
    sfilenum = -1
    openw, unit, filename, /get_lun
    save_script, unit
;    close, unit
    free_lun, unit
endif

return
bad: on_ioerror, null
if ( unit gt 0 ) then free_lun, unit
err = widget_message("IO Error has occurred writing script.")
print, "IO Error has occured writing in script."
widget_control, set_value = "Error writing script.", status_text

END
;**********************************************************

;**********************************************************
PRO save_script, unit
common shared_1
common shared_2
common shared_3
common shared_4
common shared_5
common shared_6
common shared_ao


if (aomode) then begin
 mag = 10.6
endif else begin
 mag = 1.0
endelse

flatcmd = 'flat -c 10000'

strfilter = string(boxes(curbox).filter)
strblocker = string(boxes(curbox).blocker)
;print, 'filter is ',  strfilter
;print, 'blocker is ', strblocker

if ( aomode ) then begin
  boxes(curbox).blocker=2
  if ( boxes(curbox).wheel eq 1 ) then begin
    filtercmd = 'filter ' + strfilter + ' 7'
  endif else begin
    filtercmd = 'filter 2 ' + strfilter
  endelse
endif else begin
;  widget_control, thinmenu, set_value=1
  boxes(curbox).blocker=1
  if ( boxes(curbox).wheel eq 1 ) then begin
    filtercmd = 'filter ' + strfilter + ' 11'
  endif else begin
    filtercmd = 'filter ' + strblocker + ' ' + strfilter
  endelse
endelse 

;print, filtercmd

; common to all EFS scripts
printf, unit, '#!/bin/csh -f'

;  Added in this block of code which was in place in the June 18 2011
;   version, but apparently missing starting in the Oct. 14 2011
;   version.  GWD  Dec. 11, 2012

; First write out raw parameters for later recall.
    printf, unit, format='("#",I2)', nodmode
    printf, unit, format='("#",I2)', calmode
    printf, unit, format='("#",I2)', totbox
    for i = 1, totbox do begin
        if filtermenu_index(i) eq 18 then blank_flag=1 else blank_flag=0
        printf, unit, format='("#",F12.8)', boxes(i).echphi
        printf, unit, format='("#",F12.8)', boxes(i).graphi
        printf, unit, format='("#",I2)', boxes(i).reps
        printf, unit, format='("#",A)', boxes(i).spmode
        printf, unit, format='("#",A)', boxes(i).slit
        printf, unit, format='("#",I2)', boxes(i).filter
        printf, unit, format='("#",I2)', boxes(i).blocker
        printf, unit, format='("#",I2)', boxes(i).wheel
        printf, unit, format='("#",A)', boxes(i).object
        printf, unit, format='("#",F12.8)', boxes(i).otint
        printf, unit, format='("#",A)', boxes(i).star
        printf, unit, format='("#",F12.8)', boxes(i).stint
        printf, unit, format='("#",F12.8)', boxes(i).rot
        printf, unit, format='("#",I5)', boxes(i).coadds
        printf, unit, format='("#",I5)', boxes(i).scoadds
    endfor

;  END of added block.  GWD  Dec. 11, 2012




printf, unit, '# Now begin csh commands.'
printf, unit, 'unalias m'
printf, unit, 'unalias s'
printf, unit, 'unalias w'
if ( serv eq 'nirspecsim' ) then begin
 printf, unit, 'alias s "show -s nirspecsim"'
 printf, unit, 'alias m "modify -s nirspecsim"'
 printf, unit, 'alias w "waitfor -s nirspecsim"'
endif else begin
 printf, unit, 'alias s "show -s nirspec"'
 printf, unit, 'alias m "modify -s nirspec"'
 printf, unit, 'alias w "waitfor -s nirspec"'
endelse
; First check to make sure no scripts are currently running.
printf, unit, 'm scriptdir=',scriptdir
printf, unit, 'set scriptstat = `show -terse -s nirspec script`'
printf, unit, 'if ($scriptstat == ''1'') then'
 printf, unit, '    scriptip.tk'
 printf, unit, '    exit'
printf, unit, 'endif'

; Added in lines below to prevent user from spawning a second
;  EFS script if a first one is already in progress
; GWD  Dec. 17, 2012

    printf, unit, 'm comment=''Starting script.'' '
    printf, unit, 'm script=1'
    printf, unit, 'm scriptname=',scriptname
    printf, unit, 'm firstfile=',sfilenum

; End of added block
    
printf, unit, '# NIRSPEC Script, written by NIRSPEC_EFS'
if ( aomode ) then begin
 printf, unit, '# NIRSPAO mode'
 printf, unit, '# set targwave'
; printf, unit, 'modify -s dcs2 targwave=',targwave
endif

i = 1
stri = strcompress(string(i), /remove_all)

slitoffset = 1.5/mag   
; variable to keep center beam in Nod 3 away from slit center   GWD 3/21/2012


case calmode of
  0: begin
     ; Setup only
     ; This will send mechanism moves regardless of wait4config
     printf, unit, 'set wfcc = 0'
     printf, unit, 'set wfcc = (`show -s nirspec -terse userint2`)'
     printf, unit, 'if ($wfcc == 0) then'
      printf, unit, '# First make sure Calibration unit is off.'
      printf, unit, 'm neon=0'
      printf, unit, 'm argon=0'
      printf, unit, 'm xenon=0'
      printf, unit, 'm krypton=0'
      printf, unit, 'flat off'
      printf, unit, 'm etalon=0'
      printf, unit, 'm calmpos=0   # Move out the cal mirror.'
      printf, unit, 'm calppos=0   # Move out the cal pinhole.'

;  GDopp  9/18/2014  Disabled opening up hatch by EFS
;      printf, unit, 'm calcpos=0   # Open the cal cover.'

      printf, unit, '#------------------------------'
      printf, unit, 'else'
      printf, unit, 'echo wait4config is off, stages not checked'
      printf, unit, '# wait4config is off, stages not checked'
     printf, unit, 'endif'
     printf, unit, '# Now setup instrument.'
     if filtermenu_index(i) eq 18 then blank_flag=1 else blank_flag=0
     ; Removed 2nd condition to check whether the BLANK is in  GWD 9/20/2012
     if ( lock_flag eq 0 ) then begin
      printf, unit, 'm slitpos=',boxes(i).slit
      printf, unit, 'm echlpos=',boxes(i).echphi*180.0/!pi
      printf, unit, 'm disppos=',(boxes(i).graphi*180.0/!pi+25.0)
      printf, unit, filtercmd
      printf, unit, 'wfcalm'
      printf, unit, 'wfcalp'
      printf, unit, 'wfcalc'
      printf, unit, 'wfslit'
      printf, unit, 'wfechl'
      printf, unit, 'wfdisp'
      printf, unit, 'wffilter'
     endif else begin
      printf, unit, 'echo EFS lock is on, no stage moves'
      printf, unit, '# EFS lock is on, no stage moves'
     endelse
;     printf, unit, 'scriptdone.tk'
     end
  1: begin
     ; Lamps only
     ; Checks wait4config before sending mechanisms
     printf, unit, 'set wfcc = 0'
     printf, unit, 'set wfcc = (`show -s nirspec -terse userint2`)'
     printf, unit, 'if ($wfcc == 0) then'
      printf, unit, '# First make sure Calibration unit is off.'
      printf, unit, 'm neon=0'
      printf, unit, 'm argon=0'
      printf, unit, 'm xenon=0'
      printf, unit, 'm krypton=0'
      printf, unit, 'flat off'
      printf, unit, 'm etalon=0'
      printf, unit, 'm calppos=0   # Move out the cal pinhole.'
      printf, unit, '#------------------------------'
      if filtermenu_index(i) eq 18 then blank_flag=1 else blank_flag=0
     ; Removed 2nd condition to check whether the BLANK is in  GWD 9/20/2012
      if ( lock_flag eq 0 ) then begin
       printf, unit, '# Now setup instrument.'
       printf, unit, 'm slitpos=',boxes(i).slit
       printf, unit, 'm echlpos=',boxes(i).echphi*180.0/!pi
       printf, unit, 'm disppos=',(boxes(i).graphi*180.0/!pi+25.0)
       printf, unit, filtercmd
       printf, unit, 'wfcalp'
       printf, unit, 'wfslit'
       printf, unit, 'wfechl'
       printf, unit, 'wfdisp'
       printf, unit, 'wffilter'
      endif else begin
       printf, unit, 'echo EFS lock is on, no stage moves'
       printf, unit, '# EFS lock is on, no stage moves'
      endelse
      printf, unit, 'else'
      printf, unit, 'echo wait4config is off, stages not checked'
      printf, unit, '# wait4config is off, stages not checked'
      printf, unit, 'endif'
      ;     Take calibration sequence
      printf, unit, 'm calmpos=1   # Move in the cal mirror.'
      printf, unit, 'sleep 1'
      printf, unit, 'set detval = `', flatcmd, '`'
      printf, unit, 'wfg'
      printf, unit, 'wft'
      printf, unit, 'set itime=$detval[1]'
      printf, unit, 'set coadds=$detval[2]'
      printf, unit, 'wfcalm'
      printf, unit, 'm sampmode=2'
      printf, unit, 'm itime=$itime'
      printf, unit, 'm coadds=$coadds'
      printf, unit, 'm comment= "flat lamp off" '
      printf, unit, 'goi'
      printf, unit, 'flat on'
      printf, unit, 'sleep 1'
      printf, unit, 'm comment= "lamp warming" '
      printf, unit, 'echo lamp warming'
      printf, unit, 'flatwait'
      printf, unit, 'm comment= "flat field" '
      printf, unit, 'goi'
      printf, unit, 'flat off'

; Now take arc lamp spectra            
      printf, unit, 'm neon=1'
      printf, unit, 'm comment=''Neon lamps ',stri,' '' '
      printf, unit, '#------------------------------'
      if ( boxes(i).spmode ne 'low res' ) then begin
       printf, unit, 'm argon=1'
       printf, unit, 'm xenon=1'
       printf, unit, 'm krypton=1'
       printf, unit, 'm comment=''Arc lamps ',stri,' '' '
       printf, unit, '#------------------------------'
      endif
      printf, unit, '# set up exposure values'
      ti = strtrim(string(i),1) ;
;     printf, unit, 'm object=''lamps on'' '
      if ( aomode) then begin
       printf, unit, 'm itime=5'
      endif else begin
       printf, unit, 'm itime=0.25'
       if ( boxes(i).spmode eq 'low res' ) then begin
        printf, unit, 'm itime=0.25'
       endif
      endelse
      printf, unit, 'm coadds=10'
      printf, unit, 'm sampmode=2'
      printf, unit, 'goi'

; If spectrum is low res, also take an argon lamp.  Now take arc lamp
; spectra 
      if ( boxes(i).spmode eq 'low res' ) then begin
       printf, unit, 'm neon=0'
       printf, unit, 'm argon=1'
       printf, unit, 'm comment=''Argon lamps ',stri,' '' '
       printf, unit, '#------------------------------'
       printf, unit, '# set up exposure values'
       ti = strtrim(string(i),1) ;
;      printf, unit, 'm object=''lamps on'' '
       if ( aomode) then begin
        printf, unit, 'm itime=5'
       endif else begin
        printf, unit, 'm itime=0.25'
       endelse
       printf, unit, 'm coadds=10'
       printf, unit, 'm sampmode=2'
       printf, unit, 'goi'
      endif

; Flush the detector.
       printf, unit, 'm neon=0'
       printf, unit, 'm argon=0'
       printf, unit, 'm xenon=0'
       printf, unit, 'm krypton=0'
       printf, unit, 'm calmpos=0   # Move out the cal mirror.'
;      printf, unit, 'm itime=1'
;      printf, unit, 'm coadds=2'
       printf, unit, 'm comment='' '' '
       printf, unit, 'flush'
;      printf, unit, 'm test=1'
;      printf, unit, 'wftest'
       printf, unit, 'm coadds=1'
       printf, unit, 'wfcalm'
;       printf, unit, 'scriptdone.tk'
; unknown endif below
;endif
    end

  2: begin
     ; Object only 
     if ( aomode ) then begin
      printf, unit, '# determine if we are in NGS or LGS mode'
      printf, unit, 'set aoopsmode = `show -s ao -terse aoopsmode`'
     endif
     printf, unit, 'set wfcc = 0'
     printf, unit, 'set wfcc = (`show -s nirspec -terse userint2`)'
     printf, unit, 'if ($wfcc == 0) then'
      printf, unit, '# First make sure Calibration unit is off.'
      printf, unit, 'm neon=0'
      printf, unit, 'm argon=0'
      printf, unit, 'm xenon=0'
      printf, unit, 'm krypton=0'
      printf, unit, 'flat off'
      printf, unit, 'm etalon=0'
      printf, unit, 'm calmpos=0   # Move out the cal mirror.'
      printf, unit, 'm calppos=0   # Move out the cal pinhole.'
      printf, unit, 'hatch open    # Open the hatch.'
      printf, unit, '#------------------------------'
      printf, unit, 'else'
      printf, unit, 'echo wait4config is off, stages not checked'
      printf, unit, '# wait4config is off, stages not checked'
     printf, unit, 'endif'
     for i = 1, totbox do begin
      stri = strcompress(string(i), /remove_all)
      if filtermenu_index(i) eq 18 then blank_flag=1 else blank_flag=0
      printf, unit, '# set up exposure values'
      printf, unit, 'wfg'
      printf, unit, 'wft'
      printf, unit, 'm comment=''Obs #',stri,' '' '
      ti = strtrim(string(i),1) ;
      if ( boxes(i).object ne 'default' ) then begin
       printf, unit, 'm object='' ',boxes(i).object,' '' '
      endif
      printf, unit, 'm itime=', boxes(i).otint
      printf, unit, 'm coadds=', boxes(i).coadds
      if ( boxes(i).otint ge 4 ) then begin
       printf, unit, 'm sampmode=3'
       printf, unit, 'm multispec=16'
      endif else begin
       printf, unit, 'm sampmode=2'
      endelse
      printf, unit, 'if ($wfcc == 0) then'
       if ( lock_flag eq 0 ) then begin
        printf, unit, '#------------------------------'
        printf, unit, '# Now setup instrument.'
        if blank_flag eq 0 then begin
         printf, unit, 'm slitpos=',boxes(i).slit
	 printf, unit, 'm echlpos=',boxes(i).echphi*180.0/!pi
	 printf, unit, 'm disppos=',((boxes(i).graphi+grahaang)*180.0/!pi)
        endif
        printf, unit, filtercmd
       endif else begin
        printf, unit, 'echo EFS lock is on, no stage moves.'
        printf, unit, '# EFS lock is on, no stage moves.'
       endelse
       printf, unit, 'else'
       printf, unit, 'echo wait4config is off, stages not checked.'
       printf, unit, '# wait4config is off, stages not checked.'
      printf, unit, 'endif' 
      for j = 1, boxes(i).reps do begin
       if ( nodmode eq 1 ) then begin
        printf, unit, '# Stare mode'
        printf, unit, 'wffilter'
        if blank_flag eq 0 then begin
         printf, unit, 'wfslit'
         printf, unit, 'wfechl'
         printf, unit, 'wfdisp'
        endif
        printf, unit, 'z_wfao'
        if ( aomode ) then begin
         printf, unit, 'if ( $aoopsmode == 2 ) then'
         printf, unit, '# lgsao mode, start LBWFS infinite loop'
         printf, unit, '  modify -s aolb aolbloop=1'
         printf, unit, 'endif'
        endif
        printf, unit, 'm go=1'
        printf, unit, 'wfg'
       endif
       if ( nodmode eq 2 ) then begin
        printf, unit, 'm comment=''Obs #',stri,', Exp 1.'' '
        slitlen = slitlist(boxes(i).slit).length
        if ( slitlist(boxes(i).slit).width gt $
        slitlist(boxes(i).slit).length ) then begin
         slitlen = slitlist(boxes(i).slit).width
        endif
        if ( slitlen gt 30.0/mag ) then slitlen = 30.0/mag
         printf, unit, '# nod2 mode'
         up = slitlen / 4.0
         down = 0.0 - (slitlen / 2.0)
         tu = strtrim(string(up),1) ;
         td = strtrim(string(down),1) ;
         printf, unit,'slitmove ',tu, ' 0',aomode
         printf, unit, 'centerup.tk'
         printf, unit, 'wffilter'
         printf, unit, 'wfslit'
         printf, unit, 'wfechl'
         printf, unit, 'wfdisp'
         printf, unit, 'z_wfao'
         printf, unit, 'm go=1'
         printf, unit, 'wfg'
         printf, unit,'slitmove ',td, ' 0',aomode
         printf, unit, 'centerup.tk'
         printf, unit, 'm comment=''Obs #',stri,', Exp 2.'' '
         printf, unit, 'z_wfao'
         printf, unit, 'm go=1'
         printf, unit, 'wfg'
         printf, unit,'slitmove ',tu, ' 0',aomode
        endif
        if ( nodmode eq 3 ) then begin
         printf, unit, 'm comment=''Obs #',stri, ', Exp 1.'' '
         slitlen = slitlist(boxes(i).slit).length
         if ( slitlist(boxes(i).slit).width gt $
         slitlist(boxes(i).slit).length ) then begin
          slitlen = slitlist(boxes(i).slit).width
         endif
         if ( slitlen gt 30.0/mag) then slitlen = 30.0/mag
          printf, unit, '# nod3 mode'
          up = slitlen / 2.0            ;  changed from 1.5  GWD
          down = 0.0 - (slitlen / 4.0)  ;  changed from 3.0  GWD
          tu = strtrim(string(up),1)    ;
          td = strtrim(string(down),1) ;
          off = strtrim(string(slitoffset),1)  ;  apply 1.5" offset along slit
          tdo = strtrim(string(down-slitoffset),1) ; move along slit & take out offset
          printf, unit,'slitmove ',off, ' 0',aomode
          printf, unit, 'centerup.tk'
          printf, unit, 'wffilter'
          printf, unit, 'wfslit'
          printf, unit, 'wfechl'
          printf, unit, 'wfdisp'
          printf, unit, 'z_wfao'
          
          if ( aomode ) then begin
           printf, unit, 'if ( $aoopsmode == 2 ) then'
           printf, unit, '# lgsao mode, start LBWFS infinite loop'
           printf, unit, '  modify -s aolb aolbloop=1'
           printf, unit, 'endif'
          endif
          printf, unit, 'm go=1'
          printf, unit, 'wfg'
          printf, unit,'slitmove ',td, ' 0',aomode
          printf, unit, 'centerup.tk'
          printf, unit, 'm comment=''Obs #',stri,', Exp 2.'' '
          printf, unit, 'z_wfao'
          printf, unit, 'm go=1'
          printf, unit, 'wfg'
          printf, unit,'slitmove ',tu, ' 0',aomode
          printf, unit, 'centerup.tk'
          printf, unit, 'm comment=''Obs #',stri,', Exp 3.'' '
          printf, unit, 'z_wfao'
          printf, unit, 'm go=1'
          printf, unit, 'wfg'
          printf, unit,'slitmove ',tdo, ' 0',aomode
         endif
         if ( nodmode eq 4 ) then begin
          printf, unit, 'm comment=''Obs #',stri,', Exp 1.'' '
          slitlen = slitlist(boxes(i).slit).length
          if ( slitlist(boxes(i).slit).width gt $
          slitlist(boxes(i).slit).length ) then begin
           slitlen = slitlist(boxes(i).slit).width
          endif
          if ( slitlen gt 30.0/mag ) then slitlen = 30.0/mag
           printf, unit, '# nod4 mode'
           up1 = slitlen / 8.0
           up2 = slitlen * 0.75
           down = 0.0 - (slitlen / 2.0)
           tu1 = strtrim(string(up1),1) ;
           tu2 = strtrim(string(up2),1) ;
           td = strtrim(string(down),1) ;
           printf, unit,'slitmove ',tu1, ' 0',aomode
           printf, unit, 'centerup.tk'
           printf, unit, 'wffilter'
           printf, unit, 'wfslit'
           printf, unit, 'wfechl'
           printf, unit, 'wfdisp'
           printf, unit, 'z_wfao'
           if ( aomode ) then begin
            printf, unit, 'if ( $aoopsmode == 2 ) then'
            printf, unit, '# lgsao mode, start LBWFS infinite loop'
            printf, unit, '  modify -s aolb aolbloop=1'
            printf, unit, 'endif'
           endif
           printf, unit, 'm go=1'
           printf, unit, 'wfg'
           printf, unit,'slitmove ',td, ' 0',aomode
           printf, unit, 'centerup.tk'
           printf, unit, 'm comment=''Obs #',stri,', Exp 2.'' '
           printf, unit, 'z_wfao'
           printf, unit, 'm go=1'
           printf, unit, 'wfg'
           printf, unit,'slitmove ',tu2, ' 0',aomode
           printf, unit, 'centerup.tk'
           printf, unit, 'm comment=''Obs #',stri,', Exp 3.'' '
           printf, unit, 'z_wfao'
           printf, unit, 'm go=1'
           printf, unit, 'wfg'
           printf, unit,'slitmove ',td, ' 0',aomode
           printf, unit, 'centerup.tk'
           printf, unit, 'm comment=''Obs #',stri,', Exp 4.'' '
           printf, unit, 'z_wfao'
           printf, unit, 'm go=1'
           printf, unit, 'wfg'
           printf, unit,'slitmove ',tu1, ' 0',aomode
          endif
          if ( nodmode eq 5 ) then begin
           printf, unit, 'm comment=''Obs #',stri,', Exp 1.'' '
           slitlen = slitlist(boxes(i).slit).length
           if ( slitlist(boxes(i).slit).width gt $
           slitlist(boxes(i).slit).length ) then begin
            slitlen = slitlist(boxes(i).slit).width
           endif
           if ( slitlen gt 30.0 ) then slitlen = 30.0/mag
           printf, unit, '# ABBA mode'
           up1 = slitlen / 4.0
           up2 = slitlen / 2.0
           down1 = 0.0 - (slitlen / 2.0)
           down2 = 0.0 - (slitlen / 4.0)
           tu1 = strtrim(string(up1),1) ;
           tu2 = strtrim(string(up2),1) ;
           td1 = strtrim(string(down1),1) ;
           td2 = strtrim(string(down2),1) ;
           printf, unit,'slitmove ',tu1, ' 0',aomode
           printf, unit, 'centerup.tk'
           printf, unit, 'wffilter'
           printf, unit, 'wfslit'
           printf, unit, 'wfechl'
           printf, unit, 'wfdisp'
           printf, unit, 'z_wfao'
           if ( aomode ) then begin
            printf, unit, 'if ( $aoopsmode == 2 ) then'
            printf, unit, '# lgsao mode, start LBWFS infinite loop'
            printf, unit, '  modify -s aolb aolbloop=1'
            printf, unit, 'endif'
           endif
           printf, unit, 'm go=1'
           printf, unit, 'wfg'
           printf, unit,'slitmove ',td1, ' 0',aomode
           printf, unit, 'centerup.tk'
           printf, unit, 'm comment=''Obs #',stri,', Exp 2.'' '
           printf, unit, 'z_wfao'
           printf, unit, 'm go=1'
           printf, unit, 'wfg'
           printf, unit, 'm comment=''Obs #',stri,', Exp 3.'' '
           printf, unit, 'z_wfao'
           printf, unit, 'm go=1'
           printf, unit, 'wfg'
           printf, unit,'slitmove ',tu2, ' 0',aomode
           printf, unit, 'centerup.tk'
           printf, unit, 'm comment=''Obs #',stri,', Exp 4.'' '
           printf, unit, 'z_wfao'
           printf, unit, 'm go=1'
           printf, unit, 'wfg'
           printf, unit,'slitmove ',td2, ' 0',aomode
          endif
          if ( nodmode eq 6 ) then begin
           printf, unit, 'm comment=''Obs #',stri,', Exp 1.'' '
           slitlen = slitlist(boxes(i).slit).length
           if ( slitlist(boxes(i).slit).width gt $
           slitlist(boxes(i).slit).length ) then begin
            slitlen = slitlist(boxes(i).slit).width
           endif
           if ( slitlen gt 30.0/mag ) then slitlen = 30.0/mag
           printf, unit, '# nod off mode'
           printf, unit, 'wffilter'
           printf, unit, 'wfslit'
           printf, unit, 'wfechl'
           printf, unit, 'wfdisp'
           printf, unit, 'z_wfao'
           if ( aomode ) then begin
            printf, unit, 'if ( $aoopsmode == 2 ) then'
            printf, unit, '# lgsao mode, start LBWFS infinite loop'
            printf, unit, '  modify -s aolb aolbloop=1'
            printf, unit, 'endif'
           endif
           printf, unit, 'm go=1'
           printf, unit, 'wfg'
           printf, unit, 'set node = `show -terse -s nirspec node`'
           printf, unit, 'set nodn = `show -terse -s nirspec nodn`'
           printf, unit, $
           'modify -s dcs2 raoff=$node decoff=$nodn rel2curr=t'
           printf, unit, 'm comment=''Obs #',stri,', Exp 2.'' '
           printf, unit, 'z_wfao'
           printf, unit, 'm go=1'
           printf, unit, 'wfg'
           printf, unit, 'set tnode = `echo "0.0 - $node" | bc -l`'
           printf, unit, 'set tnodn = `echo "0.0 - $nodn" | bc -l`'
           printf, unit, 'modify -s dcs2 raoff=$tnode decoff=$tnodn rel2curr=t'
           printf, unit,'sleep 3 '
          endif
          if ( nodmode eq 7 ) then begin
	   user = read_ascii(userfile)
	   userpos = user.field1
	   npos = n_elements(userpos)
           printf, unit, 'm comment= "User-defined nod"'
           printf, unit, '# User defined mode'
           printf, unit, '# ',userfile
           tu1 = strtrim(string(userpos[0]),1) ;
           printf, unit,'slitmove ',tu1, ' 0',aomode
           printf, unit, 'centerup.tk'
           printf, unit, 'wffilter'
           printf, unit, 'wfslit'
           printf, unit, 'wfechl'
           printf, unit, 'wfdisp'
           printf, unit, 'z_wfao'
           if ( aomode ) then begin
            printf, unit, 'if ( $aoopsmode == 2 ) then'
            printf, unit, '# lgsao mode, start LBWFS infinite loop'
            printf, unit, '  modify -s aolb aolbloop=1'
            printf, unit, 'endif'
           endif
           printf, unit, 'm go=1'
           printf, unit, 'wfg'
	   for k=1, npos-1 do begin
            if ( userpos[k] ne 0.0 ) then begin
             move = strtrim(string(userpos[k]),1)
	     printf, unit, 'slitmove ',move, ' 0',aomode
	     if ( k ne npos-1 ) then begin
              printf, unit, 'centerup.tk'
             endif
            endif
	    if ( k ne npos-1 ) then begin
             printf, unit, 'z_wfao'
             printf, unit, 'm go=1'
             printf, unit, 'wfg'
            endif
           endfor
          endif
         endfor
        endfor
     end

  3: begin

; Bright Object only   
;  Option added 22 March 2012  for nod pattern to be offset from science nods
;  by amount:  SLITOFFSET

     if ( aomode ) then begin
      printf, unit, '# determine if we are in NGS or LGS mode'
      printf, unit, 'set aoopsmode = `show -s ao -terse aoopsmode`'
     endif
     printf, unit, 'set wfcc = 0'
     printf, unit, 'set wfcc = (`show -s nirspec -terse userint2`)'
     printf, unit, 'if ($wfcc == 0) then'
      printf, unit, '# First make sure Calibration unit is off.'
      printf, unit, 'm neon=0'
      printf, unit, 'm argon=0'
      printf, unit, 'm xenon=0'
      printf, unit, 'm krypton=0'
      printf, unit, 'flat off'
      printf, unit, 'm etalon=0'
      printf, unit, 'm calmpos=0   # Move out the cal mirror.'
      printf, unit, 'm calppos=0   # Move out the cal pinhole.'
      printf, unit, 'hatch open    # Open the hatch.'
      printf, unit, '#------------------------------'
      printf, unit, 'else'
      printf, unit, 'echo wait4config is off, stages not checked'
      printf, unit, '# wait4config is off, stages not checked'
     printf, unit, 'endif'
     for i = 1, totbox do begin
      stri = strcompress(string(i), /remove_all)
      if filtermenu_index(i) eq 18 then blank_flag=1 else blank_flag=0
      printf, unit, '# set up exposure values'
      printf, unit, 'wfg'
      printf, unit, 'wft'
      printf, unit, 'm comment=''Obs #',stri,' '' '
      ti = strtrim(string(i),1) ;
      if ( boxes(i).object ne 'default' ) then begin
       printf, unit, 'm object='' ',boxes(i).object,' '' '
      endif
      printf, unit, 'm itime=', boxes(i).otint
      printf, unit, 'm coadds=', boxes(i).coadds
      if ( boxes(i).otint ge 4 ) then begin
       printf, unit, 'm sampmode=3'
       printf, unit, 'm multispec=16'
      endif else begin
       printf, unit, 'm sampmode=2'
      endelse
      printf, unit, 'if ($wfcc == 0) then'
       if ( lock_flag eq 0 ) then begin
        printf, unit, '#------------------------------'
        printf, unit, '# Now setup instrument.'
        if blank_flag eq 0 then begin
         printf, unit, 'm slitpos=',boxes(i).slit
	 printf, unit, 'm echlpos=',boxes(i).echphi*180.0/!pi
	 printf, unit, 'm disppos=',((boxes(i).graphi+grahaang)*180.0/!pi)
        endif
        printf, unit, filtercmd
       endif else begin
        printf, unit, 'echo EFS lock is on, no stage moves.'
        printf, unit, '# EFS lock is on, no stage moves.'
       endelse
       printf, unit, 'else'
       printf, unit, 'echo wait4config is off, stages not checked.'
       printf, unit, '# wait4config is off, stages not checked.'
      printf, unit, 'endif' 
      for j = 1, boxes(i).reps do begin
       if ( nodmode eq 1 ) then begin
        printf, unit, '# Stare mode'
        printf, unit, 'wffilter'
        if blank_flag eq 0 then begin
         printf, unit, 'wfslit'
         printf, unit, 'wfechl'
         printf, unit, 'wfdisp'
        endif
        printf, unit, 'z_wfao'
        if ( aomode ) then begin
         printf, unit, 'if ( $aoopsmode == 2 ) then'
         printf, unit, '# lgsao mode, start LBWFS infinite loop'
         printf, unit, '  modify -s aolb aolbloop=1'
         printf, unit, 'endif'
        endif
        printf, unit, 'm go=1'
        printf, unit, 'wfg'
       endif
       if ( nodmode eq 2 ) then begin
        printf, unit, 'm comment=''Obs #',stri,', Exp 1.'' '
        slitlen = slitlist(boxes(i).slit).length
        if ( slitlist(boxes(i).slit).width gt $
        slitlist(boxes(i).slit).length ) then begin
         slitlen = slitlist(boxes(i).slit).width
        endif
        if ( slitlen gt 30.0/mag ) then slitlen = 30.0/mag
         printf, unit, '# nod2 mode'
         up = slitlen / 4.0
         down = 0.0 - (slitlen / 2.0)
         tuo = strtrim(string(up+slitoffset),1) ; add 1.5" offset along slit
         tuo2 = strtrim(string(up-slitoffset),1) ; sub. 1.5" offset along slit
         td = strtrim(string(down),1) ;
         printf, unit,'slitmove ',tuo, ' 0',aomode
         printf, unit, 'centerup.tk'
         printf, unit, 'wffilter'
         printf, unit, 'wfslit'
         printf, unit, 'wfechl'
         printf, unit, 'wfdisp'
         printf, unit, 'z_wfao'
         printf, unit, 'm go=1'
         printf, unit, 'wfg'
         printf, unit,'slitmove ',td, ' 0',aomode
         printf, unit, 'centerup.tk'
         printf, unit, 'm comment=''Obs #',stri,', Exp 2.'' '
         printf, unit, 'z_wfao'
         printf, unit, 'm go=1'
         printf, unit, 'wfg'
         printf, unit,'slitmove ',tuo2,' 0',aomode
        endif
        if ( nodmode eq 3 ) then begin
         printf, unit, 'm comment=''Obs #',stri, ', Exp 1.'' '
         slitlen = slitlist(boxes(i).slit).length
         if ( slitlist(boxes(i).slit).width gt $
         slitlist(boxes(i).slit).length ) then begin
          slitlen = slitlist(boxes(i).slit).width
         endif
         if ( slitlen gt 30.0/mag) then slitlen = 30.0/mag
          printf, unit, '# nod3 mode'
          up = slitlen / 2.0            ;  changed from 1.5  GWD
          down = 0.0 - (slitlen / 4.0)  ;  changed from 3.0  GWD
          tu = strtrim(string(up),1)    ;
          td = strtrim(string(down),1) ;
          off = strtrim(string(-slitoffset),1)  ;  apply 1.5" offset along slit
          tdo = strtrim(string(down+slitoffset),1) ; move along slit & take out offset
          printf, unit,'slitmove ',off, ' 0',aomode
          printf, unit, 'centerup.tk'
          printf, unit, 'wffilter'
          printf, unit, 'wfslit'
          printf, unit, 'wfechl'
          printf, unit, 'wfdisp'
          printf, unit, 'z_wfao'
          
          if ( aomode ) then begin
           printf, unit, 'if ( $aoopsmode == 2 ) then'
           printf, unit, '# lgsao mode, start LBWFS infinite loop'
           printf, unit, '  modify -s aolb aolbloop=1'
           printf, unit, 'endif'
          endif
          printf, unit, 'm go=1'
          printf, unit, 'wfg'
          printf, unit,'slitmove ',td, ' 0',aomode
          printf, unit, 'centerup.tk'
          printf, unit, 'm comment=''Obs #',stri,', Exp 2.'' '
          printf, unit, 'z_wfao'
          printf, unit, 'm go=1'
          printf, unit, 'wfg'
          printf, unit,'slitmove ',tu, ' 0',aomode
          printf, unit, 'centerup.tk'
          printf, unit, 'm comment=''Obs #',stri,', Exp 3.'' '
          printf, unit, 'z_wfao'
          printf, unit, 'm go=1'
          printf, unit, 'wfg'
          printf, unit,'slitmove ',tdo, ' 0',aomode
         endif
         if ( nodmode eq 4 ) then begin
          printf, unit, 'm comment=''Obs #',stri,', Exp 1.'' '
          slitlen = slitlist(boxes(i).slit).length
          if ( slitlist(boxes(i).slit).width gt $
          slitlist(boxes(i).slit).length ) then begin
           slitlen = slitlist(boxes(i).slit).width
          endif
          if ( slitlen gt 30.0/mag ) then slitlen = 30.0/mag
           printf, unit, '# nod4 mode'
           up1 = slitlen / 8.0
           up2 = slitlen * 0.75
           down = 0.0 - (slitlen / 2.0)
           tu1 = strtrim(string(up1),1) ;
           tu2 = strtrim(string(up2),1) ;
           td = strtrim(string(down),1) ;
           printf, unit,'slitmove ',tu1, ' 0',aomode
           printf, unit, 'centerup.tk'
           printf, unit, 'wffilter'
           printf, unit, 'wfslit'
           printf, unit, 'wfechl'
           printf, unit, 'wfdisp'
           printf, unit, 'z_wfao'
           if ( aomode ) then begin
            printf, unit, 'if ( $aoopsmode == 2 ) then'
            printf, unit, '# lgsao mode, start LBWFS infinite loop'
            printf, unit, '  modify -s aolb aolbloop=1'
            printf, unit, 'endif'
           endif
           printf, unit, 'm go=1'
           printf, unit, 'wfg'
           printf, unit,'slitmove ',td, ' 0',aomode
           printf, unit, 'centerup.tk'
           printf, unit, 'm comment=''Obs #',stri,', Exp 2.'' '
           printf, unit, 'z_wfao'
           printf, unit, 'm go=1'
           printf, unit, 'wfg'
           printf, unit,'slitmove ',tu2, ' 0',aomode
           printf, unit, 'centerup.tk'
           printf, unit, 'm comment=''Obs #',stri,', Exp 3.'' '
           printf, unit, 'z_wfao'
           printf, unit, 'm go=1'
           printf, unit, 'wfg'
           printf, unit,'slitmove ',td, ' 0',aomode
           printf, unit, 'centerup.tk'
           printf, unit, 'm comment=''Obs #',stri,', Exp 4.'' '
           printf, unit, 'z_wfao'
           printf, unit, 'm go=1'
           printf, unit, 'wfg'
           printf, unit,'slitmove ',tu1, ' 0',aomode
          endif
          if ( nodmode eq 5 ) then begin
           printf, unit, 'm comment=''Obs #',stri,', Exp 1.'' '
           slitlen = slitlist(boxes(i).slit).length
           if ( slitlist(boxes(i).slit).width gt $
           slitlist(boxes(i).slit).length ) then begin
            slitlen = slitlist(boxes(i).slit).width
           endif
           if ( slitlen gt 30.0 ) then slitlen = 30.0/mag
           printf, unit, '# ABBA mode'
           up1 = slitlen / 4.0
           up2 = slitlen / 2.0
           down1 = 0.0 - (slitlen / 2.0)
           down2 = 0.0 - (slitlen / 4.0)
           tu1 = strtrim(string(up1+slitoffset),1) ;
           tu2 = strtrim(string(up2),1) ;
           td1 = strtrim(string(down1),1) ;
           td2 = strtrim(string(down2-slitoffset),1) ;
           printf, unit,'slitmove ',tu1, ' 0',aomode
           printf, unit, 'centerup.tk'
           printf, unit, 'wffilter'
           printf, unit, 'wfslit'
           printf, unit, 'wfechl'
           printf, unit, 'wfdisp'
           printf, unit, 'z_wfao'
           if ( aomode ) then begin
            printf, unit, 'if ( $aoopsmode == 2 ) then'
            printf, unit, '# lgsao mode, start LBWFS infinite loop'
            printf, unit, '  modify -s aolb aolbloop=1'
            printf, unit, 'endif'
           endif
           printf, unit, 'm go=1'
           printf, unit, 'wfg'
           printf, unit,'slitmove ',td1, ' 0',aomode
           printf, unit, 'centerup.tk'
           printf, unit, 'm comment=''Obs #',stri,', Exp 2.'' '
           printf, unit, 'z_wfao'
           printf, unit, 'm go=1'
           printf, unit, 'wfg'
           printf, unit, 'm comment=''Obs #',stri,', Exp 3.'' '
           printf, unit, 'z_wfao'
           printf, unit, 'm go=1'
           printf, unit, 'wfg'
           printf, unit,'slitmove ',tu2, ' 0',aomode
           printf, unit, 'centerup.tk'
           printf, unit, 'm comment=''Obs #',stri,', Exp 4.'' '
           printf, unit, 'z_wfao'
           printf, unit, 'm go=1'
           printf, unit, 'wfg'
           printf, unit,'slitmove ',td2, ' 0',aomode
          endif
          if ( nodmode eq 6 ) then begin
           printf, unit, 'm comment=''Obs #',stri,', Exp 1.'' '
           slitlen = slitlist(boxes(i).slit).length
           if ( slitlist(boxes(i).slit).width gt $
           slitlist(boxes(i).slit).length ) then begin
            slitlen = slitlist(boxes(i).slit).width
           endif
           if ( slitlen gt 30.0/mag ) then slitlen = 30.0/mag
           printf, unit, '# nod off mode'
           printf, unit, 'wffilter'
           printf, unit, 'wfslit'
           printf, unit, 'wfechl'
           printf, unit, 'wfdisp'
           printf, unit, 'z_wfao'
           if ( aomode ) then begin
            printf, unit, 'if ( $aoopsmode == 2 ) then'
            printf, unit, '# lgsao mode, start LBWFS infinite loop'
            printf, unit, '  modify -s aolb aolbloop=1'
            printf, unit, 'endif'
           endif
           printf, unit, 'm go=1'
           printf, unit, 'wfg'
           printf, unit, 'set node = `show -terse -s nirspec node`'
           printf, unit, 'set nodn = `show -terse -s nirspec nodn`'
           printf, unit, $
           'modify -s dcs2 raoff=$node decoff=$nodn rel2curr=t'
           printf, unit, 'm comment=''Obs #',stri,', Exp 2.'' '
           printf, unit, 'z_wfao'
           printf, unit, 'm go=1'
           printf, unit, 'wfg'
           printf, unit, 'set tnode = `echo "0.0 - $node" | bc -l`'
           printf, unit, 'set tnodn = `echo "0.0 - $nodn" | bc -l`'
           printf, unit, 'modify -s dcs2 raoff=$tnode decoff=$tnodn rel2curr=t'
           printf, unit,'sleep 3 '
          endif
          if ( nodmode eq 7 ) then begin
	   user = read_ascii(userfile)
	   userpos = user.field1
	   npos = n_elements(userpos)
           printf, unit, 'm comment= "User-defined nod"'
           printf, unit, '# User defined mode'
           printf, unit, '# ',userfile
           tu1 = strtrim(string(userpos[0]),1) ;
           printf, unit,'slitmove ',tu1, ' 0',aomode
           printf, unit, 'centerup.tk'
           printf, unit, 'wffilter'
           printf, unit, 'wfslit'
           printf, unit, 'wfechl'
           printf, unit, 'wfdisp'
           printf, unit, 'z_wfao'
           if ( aomode ) then begin
            printf, unit, 'if ( $aoopsmode == 2 ) then'
            printf, unit, '# lgsao mode, start LBWFS infinite loop'
            printf, unit, '  modify -s aolb aolbloop=1'
            printf, unit, 'endif'
           endif
           printf, unit, 'm go=1'
           printf, unit, 'wfg'
	   for k=1, npos-1 do begin
            if ( userpos[k] ne 0.0 ) then begin
             move = strtrim(string(userpos[k]),1)
	     printf, unit, 'slitmove ',move, ' 0',aomode
	     if ( k ne npos-1 ) then begin
              printf, unit, 'centerup.tk'
             endif
            endif
	    if ( k ne npos-1 ) then begin
             printf, unit, 'z_wfao'
             printf, unit, 'm go=1'
             printf, unit, 'wfg'
            endif
           endfor
          endif
         endfor
        endfor
     end

else: begin
     ; Unknown or deprecated modes...exit
      end
endcase
printf, unit, 'm script=0'
printf, unit, 'm scriptname=NULL'
printf, unit, 'm comment=''Script is finished.'' '
printf, unit, 'scriptdone.tk'

end

;**********************************************************
PRO new_box_event, Event
common shared_1
common shared_2
common shared_3
common shared_4
common shared_5
common shared_6
common shared_ao

; increase total number of boxes
if (totbox lt (maxbox-1)) then begin

; if statement for arrow keys for what to pass to calc_pixel
if in_draw then begin  ; if cursor is in draw window
    pix_x=csr_x
    pix_y=csr_y
endif else begin
    pix_x=box_x(curbox)
    pix_y=box_y(curbox)
endelse

    totbox = totbox + 1
    for i = totbox, curbox+1, -1 do begin
        boxes(i) = boxes(i-1)
	filtermenu_index(i)=filtermenu_index(i-1)
	box_x(i)=box_x(i-1)
	box_y(i)=box_y(i-1)
    endfor
    ; increment current box to be new box
    curbox = curbox + 1
    ; fill initial values
    filtermenu_index(curbox)=current_fil
    boxes(curbox).filter = filterarr(current_fil).position
    boxes(curbox).wheel = filterarr(current_fil).wheel
    boxes(curbox).slit=current_slit
    boxes(curbox).spmode=sptype
if sptype eq 'imaging'then begin
; Changed temporarily for first light to 179.99 degrees
    boxes(curbox).echphi = roundphi(179.7*!dtor)
;    boxes(curbox).echphi=roundphi(!pi)
    boxes(curbox).graphi=0.0
    box_x(curbox)=box_x(curbox-1)
    box_y(curbox)=box_y(curbox-1)
endif else begin
    boxes(curbox).echphi = roundphi(echtheta)
    boxes(curbox).graphi = roundphi(graphi(499))
    ; draw new box
    if mask_flag then erase_bad
    find_corners

; enter center of box in data coords
tmp=[(boxes(curbox).x(0)+boxes(curbox).x(2))/2.0, $
	(boxes(curbox).y(0)+boxes(curbox).y(2))/2.0]
temp=convert_coord(tmp(0), tmp(1), /data, /to_device)
box_x(curbox)=temp(0)
box_y(curbox)=temp(1)
calc_pixel, pix_x, pix_y

endelse

    update_buttons
    draw_echellogram
    draw_boxes

endif

END
;**********************************************************

;**********************************************************
PRO del_box_event, Event
common shared_1
common shared_2
common shared_3
common shared_4
common shared_5
common shared_6
common shared_ao

; decrement total number of boxes
if (totbox gt 1 and curbox ne 0) then begin

; if statement for arrow keys for what to pass to calc_pixel
if in_draw then begin
    pix_x=csr_x
    pix_y=csr_y
endif else begin
    pix_x=box_x(curbox)
    pix_y=box_y(curbox)
endelse
    if curbox eq totbox then begin
        curbox=curbox-1
    endif else begin
        for i = curbox, totbox-1 do begin
            boxes(i) = boxes(i+1)
	    filtermenu_index(i)=filtermenu_index(i+1)
	    box_x(i)=box_x(i+1)
	    box_y(i)=box_y(i+1)
        endfor
    endelse

    totbox = totbox - 1

   slitlen = slitlist(boxes(curbox).slit).length
    update_buttons
    calc_pixel, pix_x, pix_y
    choose_standard_filter
    setup_ech_filter
    if mask_flag then erase_bad
    find_corners
    draw_echellogram
    draw_boxes

endif

END
;**********************************************************

;**********************************************************
PRO next_box_event, Event
common shared_1
common shared_2
common shared_3
common shared_4
common shared_5
common shared_6
common shared_ao

if (curbox lt totbox) then begin

; if statement for arrow keys for what to pass to calc_pixel
if in_draw then begin
    pix_x=csr_x
    pix_y=csr_y
endif else begin
    pix_x=box_x(curbox)
    pix_y=box_y(curbox)
endelse

    curbox = curbox + 1
    slitlen = slitlist(boxes(curbox).slit).length
    update_buttons
    calc_pixel, pix_x, pix_y
    choose_standard_filter
    setup_ech_filter
    if mask_flag then erase_bad
    find_corners
    draw_echellogram
    draw_boxes

endif

END
;**********************************************************

;**********************************************************
PRO pre_box_event, Event
common shared_1
common shared_2
common shared_3
common shared_4
common shared_5
common shared_6
common shared_ao

if (curbox gt 1) then begin

; if statement for arrow keys for what to pass to calc_pixel
if in_draw then begin
    pix_x=csr_x
    pix_y=csr_y
endif else begin
    pix_x=box_x(curbox)
    pix_y=box_y(curbox)
endelse

    curbox = curbox - 1
    slitlen = slitlist(boxes(curbox).slit).length
    update_buttons
    calc_pixel, pix_x, pix_y
    choose_standard_filter
    setup_ech_filter
    if mask_flag then erase_bad
    find_corners
    draw_echellogram
    draw_boxes

endif

END
;**********************************************************

;**********************************************************
pro curbox_button_event, event
; displays information about the current box.
common shared_1
common shared_2
common shared_3
common shared_4
common shared_5
common shared_6
common shared_ao


if calmode eq 0 or calmode eq 1 then begin 
    objtmp='N/A'
    tinttmp='N/A'
    coaddstmp='N/A'
endif else begin
    objtmp=boxes(curbox).object
    tinttmp=strcompress(boxes(curbox).otint, /remove_all)
    coaddstmp=strcompress(boxes(curbox).coadds, /remove_all)
endelse

if calmode eq 4 or calmode eq 5 then begin
    startmp=boxes(curbox).star
    stinttmp=strcompress(boxes(curbox).stint, /remove_all)
    scoaddstmp=strcompress(boxes(curbox).scoadds, /remove_all)
endif else begin
    startmp='N/A'
    stinttmp='N/A'
    scoaddstmp='N/A'
endelse

if boxes(curbox).spmode eq 'low res' then echtmp='179.7' else $
  echtmp=strcompress(round2((boxes(curbox).echphi)*180/!pi), /remove_all)

title='INFORMATION ABOUT CURRENT BOX: '+strcompress(curbox)
objtx    ='   OBJECT NAME:            '+objtmp
tinttx   ='   INTEGRATION TIME:       '+tinttmp
coaddstx ='   OBJECT COADDS:          '+coaddstmp

startx   ='   STAR NAME:              '+startmp
stinttx  ='   INTEGRATION TIME:       '+stinttmp
scoaddstx='   STAR COADDS:            '+scoaddstmp

spmodetx ='   SPECTRAL MODE:          '+boxes(curbox).spmode
repstx   ='   NUMBER OF REPS:        '+strcompress(boxes(curbox).reps)
filtx    ='   FILTER:                 '+filterarr[filtermenu_index(curbox)].name
blocktx  ='   BLOCKING:               '+thinname[boxes(curbox).blocker]
slittx   ='   SLIT:                   '+slitlist[boxes(curbox).slit].name
echtx    ='   ECHELLE ANGLE           '+echtmp 
gratx    ='   CROSS DISPERSOR:       '+ $
  strcompress(round2((boxes(curbox).graphi+grahaang)*180/!pi))

text=[title, '', objtx, tinttx, coaddstx, '', startx, stinttx, scoaddstx, $
      '', spmodetx, repstx, filtx, blocktx, slittx, echtx, gratx]

junk=dialog_message(text, /information, dialog_parent=event.id, $
                    title='CURRENT BOX: '+strcompress(curbox))

end

;**********************************************************

;**********************************************************
PRO boxnum_event, Event
; Handles entries to the boxnum box.
common shared_1
common shared_2
common shared_3
common shared_4
common shared_5
common shared_6
common shared_ao

widget_control, Event.id, get_value= ctmp
i = strlen(ctmp(0))
if ( i gt 0 ) then begin
    reads, ctmp, i
    i = fix(i)
    if ( (i gt 0) and (i le totbox) ) then begin
        curbox = i
	slitlen = slitlist(boxes(curbox).slit).length
	update_buttons
	calc_pixel, box_x(curbox), box_y(curbox)
        choose_standard_filter
        setup_ech_filter
        if mask_flag then erase_bad
        find_corners
        draw_echellogram
        draw_boxes


    endif
endif

END
;**********************************************************

;**********************************************************
PRO itime_event, Event
common shared_1
common shared_2
common shared_3
common shared_4
common shared_5
common shared_6
common shared_ao

f = 1.0
min = 0.25
;ctmp = 3.5
widget_control, Event.id, get_value=ctmp
i = strlen(ctmp(0))
if ( i gt 0 ) then begin
    reads, ctmp, f
    if ( f lt min ) then begin
      f = min
    endif
    boxes(curbox).otint(0) = f
endif
END
;**********************************************************

;**********************************************************
PRO grabox_event, Event
common shared_1
common shared_2
common shared_3
common shared_4
common shared_5
common shared_6
common shared_ao

f = 1.0
;ctmp = 3.5
atmp = 'hi'
widget_control, Event.id, get_value=atmp
i = strlen(atmp(0))
if ( i gt 0 ) then begin
    reads, atmp, f
    boxes(curbox).graphi = roundphi(f*!pi/180.0 - grahaang)
    erase_curbox                ; Erase current box
    find_corners
    erase_curbox                ; Draw current box
    if mask_flag then begin
        erase_bad
        find_bad
        erase_bad
    endif

endif

END
;**********************************************************


;**********************************************************
PRO echbox_event, Event
common shared_1
common shared_2
common shared_3
common shared_4
common shared_5
common shared_6
common shared_ao

f = 1.0
;ctmp = 3.5
ctmp = 'hi'
widget_control, Event.id, get_value=ctmp
i = strlen(ctmp(0))
if ( i gt 0 ) then begin
    reads, ctmp, f
    boxes(curbox).echphi = roundphi(f*!pi/180.0)
    erase_curbox                ; Erase current box
    find_corners
    erase_curbox                ; Draw current box
    if mask_flag then begin
        erase_bad
        find_bad
        erase_bad
    endif

endif

END
;**********************************************************

;**********************************************************
PRO stime_event, Event
common shared_1
common shared_2
common shared_3
common shared_4
common shared_5
common shared_6
common shared_ao

ctmp = 3.5
widget_control, Event.id, get_value=ctmp
i = strlen(ctmp(0))
if ( i gt 0 ) then begin
    reads, ctmp, i
    boxes(curbox).stint(0) = i
endif

END
;**********************************************************

;**********************************************************
PRO coadds_event, Event
common shared_1
common shared_2
common shared_3
common shared_4
common shared_5
common shared_6
common shared_ao

ctmp = 1
widget_control, Event.id, get_value=ctmp
i = strlen(ctmp(0))
if ( i gt 0 ) then begin
    reads, ctmp, i
    boxes(curbox).coadds(0) = fix(i)
endif

END
;**********************************************************

;**********************************************************
PRO scoadds_event, Event
common shared_1
common shared_2
common shared_3
common shared_4
common shared_5
common shared_6
common shared_ao

ctmp = 1
widget_control, Event.id, get_value=ctmp
i = strlen(ctmp(0))
if ( i gt 0 ) then begin
    reads, ctmp, i
    boxes(curbox).scoadds(0) = fix(i)
endif

END
;**********************************************************

;**********************************************************
PRO obj_name_event, Event
common shared_1
common shared_2
common shared_3
common shared_4
common shared_5
common shared_6
common shared_ao

name = 'hi there'
widget_control, Event.id, get_value=name
boxes(curbox).object(0) = name

END
;**********************************************************

;**********************************************************
PRO star_name_event, Event
common shared_1
common shared_2
common shared_3
common shared_4
common shared_5
common shared_6
common shared_ao

tname = 'hi there'
widget_control, Event.id, get_value=tname
boxes(curbox).star(0) = tname

END
;**********************************************************

;**********************************************************
PRO update_buttons
common shared_1
common shared_2
common shared_3
common shared_4
common shared_5
common shared_6
common shared_ao
common shared_go

; Set the text boxes.
widget_control, set_value = boxes(curbox).object, obj_name
widget_control, set_value = boxes(curbox).star, star_name
widget_control, set_value = boxes(curbox).otint, itime
widget_control, set_value = boxes(curbox).stint, stime
widget_control, set_value = boxes(curbox).coadds, coadds
if not (aomode) then begin
 widget_control, thinmenu, set_value = boxes(curbox).blocker
endif
widget_control, set_value = strcompress(curbox), curbox_field
widget_control, set_value = 'Ready to set up instrument.', status_text

; Set the spectroscopy mode.
if ( boxes(curbox).spmode eq 'high res' ) then begin
    widget_control, set_droplist_select = 0, specmode
endif
if ( boxes(curbox).spmode eq 'low res' ) then begin
    widget_control, set_droplist_select = 1, specmode
endif
if ( boxes(curbox).spmode eq 'imaging' ) then begin
    widget_control, set_droplist_select = 2, specmode
endif
sptype=boxes(curbox).spmode

; Set the filter
;for i = 0, (numfil-1) do begin
;    if ( (boxes(curbox).filter eq filterarr(i).position) and $
;         (boxes(curbox).wheel eq filterarr(i).wheel) ) then begin
;        widget_control, set_droplist_select = i, filtermenu
;	current_fil=i
;    endif
;endfor
current_fil=filtermenu_index(curbox)
widget_control, filtermenu, set_droplist_select=current_fil

; disable box buttons as necessary
;if (calmode eq 0) or (calmode eq 3) or (calmode eq 4) then begin
;  new=0

if (totbox eq maxbox) or (calmode eq 0) or (calmode eq 1) or (calmode eq 4) $
  or (calmode eq 5) then new=0 else new=1
if curbox eq 1 then pre=0 else pre=1
if curbox eq totbox then next=0 else next=1
if totbox eq 1 then del=0 else del=1

widget_control, new_box, sensitive=new
widget_control, pre_box, sensitive=pre
widget_control, next_box, sensitive=next
widget_control, del_box, sensitive=del

; disable blocking menu if filter in 1st wheel
if not (aomode) then begin
 if filterarr(filtermenu_index(curbox)).wheel eq 1 then begin
;     widget_control, thinmenu, set_value=4
     widget_control, thinmenu, sensitive=0
 endif else begin
     widget_control, thinmenu, sensitive=1
 endelse
endif

; disable filter plotting for PaBeta, HeI, and M
pa_beta=where(filterarr.name eq 'PA-BETA(1.28)')
hei=where(filterarr.name eq 'HeI(1.08)')
blk=where(filterarr.name eq 'BLANK')
if (filtermenu_index(curbox) eq pa_beta(0)) or $
  (filtermenu_index(curbox) eq hei(0)) or $
  (filtermenu_index(curbox) eq blk(0)) then begin
    widget_control, plot_filter_button, sensitive=0
endif else begin
    widget_control, plot_filter_button, sensitive=1
endelse


widget_control, set_droplist_select = ((boxes(curbox).reps)-1), reps
widget_control, set_droplist_select = calmode, full
widget_control, set_droplist_select = (nodmode-1), nods

; Set the slit width
widget_control, set_droplist_select = boxes(curbox).slit, slitmenu
current_slit=boxes(curbox).slit

END
;**********************************************************
pro cfg_scr_buttons_event, event
common shared_1
common shared_2
common shared_3
common shared_4
common shared_5
common shared_6
common shared_ao

widget_control, event.top, get_uval=cfg_uval

widget_control, cfg_uval.dirbox, get_value=cfg_dir
widget_control, cfg_uval.namebox, get_value=cfg_name

case event.value of
	'OK': begin 
		num_scr=fix(cfg_name[0])
		scriptdir=cfg_dir[0]
		widget_control, event.top, /destroy
	end
	'APPLY': begin
		num_scr=fix(cfg_name[0])
		scriptdir=cfg_dir[0]
	end
	'CANCEL': widget_control, event.top, /destroy
endcase

end
;**********************************************************

;**********************************************************
pro config_script
common shared_1
common shared_2
common shared_3
common shared_4
common shared_5
common shared_6
common shared_ao

scr_name='Script Name:     '+scr_name_base
scr_num=string(format='(I4.4)',num_scr)
cfg_scr_base=widget_base(/col, title='Configure Script File')
cfg_scr_dir_box=cw_field(cfg_scr_base, value=scriptdir, font=font, $
	title='Script Directory:')
cfg_name_base=widget_base(cfg_scr_base, /row)
cfg_scr_name_box=cw_field(cfg_name_base, value=scr_num, font=font, $ 
	title=scr_name, xs=4)
cfg_scr_name_ext=widget_label(cfg_name_base, value='.csh', font=font)
cfg_scr_buttons=cw_bgroup(cfg_scr_base, ['OK', 'APPLY', 'CANCEL'], /row, $
	/return_name, font=font)

uval={dirbox:cfg_scr_dir_box, namebox:cfg_scr_name_box}

widget_control, cfg_scr_base, /realize, set_uval=uval

xmanager, 'cfg_scr_buttons', cfg_scr_buttons, /just_reg, /no_block
xmanager

END
;**********************************************************



;**********************************************************
PRO SAVE_CONFIG
common shared_1
common shared_2
common shared_3
common shared_5
common shared_4
common shared_6
common shared_ao


END
;**********************************************************

;**********************************************************
PRO OPEN_CONFIG
common shared_1
common shared_2
common shared_3
common shared_4
common shared_5
common shared_6
common shared_ao

on_ioerror, bad
unit = 0
tmpchar = 'hi'
new_path=''
filename = dialog_pickfile(filter='*', get_path=new_path, path=working_dir)
j = strlen(filename)
if ( j gt 0 ) then begin
    working_dir=new_path
    openr, unit, filename, /get_lun
    readf, unit, tmpchar
; First read out raw parameters.
    readf, unit, format='(A1,I2)', tmpchar, nodmode
    readf, unit, format='(A1,I2)', tmpchar, calmode
    readf, unit, format='(A1,I2)', tmpchar, totbox
    for i = 1, totbox do begin
        readf, unit, format='(A)', tmpchar
        val = strmid(tmpchar,1,20)
        boxes(i).echphi = float( strtrim(val,1))
        readf, unit, format='(A)', tmpchar
        val = strmid(tmpchar,1,20)
        boxes(i).graphi = float( strtrim(val,1))
        readf, unit, format='(A)', tmpchar
        val = strmid(tmpchar,1,20)
        boxes(i).reps = fix( strtrim(val,1))
        readf, unit, format='(A)', tmpchar
        val = strmid(tmpchar,1,30)
        boxes(i).spmode = strtrim(val)
        readf, unit, format='(A)', tmpchar
        val = strmid(tmpchar,1,30)
        boxes(i).slit = strtrim(val)
        readf, unit, format='(A)', tmpchar
        val = strmid(tmpchar,1,30)
        boxes(i).filter = fix(strtrim(val))
        readf, unit, format='(A)', tmpchar
        val = strmid(tmpchar,1,20)
        boxes(i).blocker = fix( strtrim(val,1))
        readf, unit, format='(A)', tmpchar
        val = strmid(tmpchar,1,20)
        boxes(i).wheel = fix( strtrim(val,1))
        readf, unit, format='(A)', tmpchar
        val = strmid(tmpchar,1,30)
        boxes(i).object = strtrim(val)
        readf, unit, format='(A)', tmpchar
        val = strmid(tmpchar,1,20)
        boxes(i).otint = float( strtrim(val,1))
        readf, unit, format='(A)', tmpchar
        val = strmid(tmpchar,1,30)
        boxes(i).star = strtrim(val)
        readf, unit, format='(A)', tmpchar
        val = strmid(tmpchar,1,20)
        boxes(i).stint = float( strtrim(val,1))
        readf, unit, format='(A)', tmpchar
        val = strmid(tmpchar,1,20)
        boxes(i).rot = float( strtrim(val,1))
        readf, unit, format='(A)', tmpchar
        val = strmid(tmpchar,1,20)
        boxes(i).coadds = fix(strtrim(val,1))
        readf, unit, format='(A)', tmpchar
        val = strmid(tmpchar,1,20)
        boxes(i).scoadds = fix(strtrim(val,1))
        curbox = i
        choose_standard_filter
        find_corners
    endfor
;    close, unit
    free_lun, unit
    curbox = totbox
    slitlen = slitlist(boxes(curbox).slit).length
    fill_curbox_arrays
    update_buttons
    choose_standard_filter
    setup_ech_filter
    find_corners
    draw_echellogram
    draw_boxes
    case calmode of 
        0:begin                ; Setup Only
            widget_control, obj_name, sensitive=0
            widget_control, star_name, sensitive=0
            widget_control, itime, sensitive=0
            widget_control, stime, sensitive=0
            widget_control, coadds, sensitive=0
            widget_control, scoadds, sensitive=0
;            widget_control, boxnum, sensitive=0
            widget_control, box_button_base, sensitive=0
            widget_control, i_boxbase, sensitive=0
            curbox=1
            totbox=1
;            widget_control, boxnum, set_value=1
            widget_control, curbox_field, set_value=strcompress(1)
            update_buttons
            calc_pixel, box_x(curbox), box_y(curbox)
            choose_standard_filter
            setup_ech_filter
            if mask_flag then erase_bad
            find_corners
            draw_echellogram
            draw_boxes
        end
        1:begin                 ; Lamps Only
            widget_control, obj_name, sensitive=0
            widget_control, star_name, sensitive=0
            widget_control, itime, sensitive=0
            widget_control, stime, sensitive=0
            widget_control, coadds, sensitive=0
            widget_control, scoadds, sensitive=0
;           widget_control, boxnum, sensitive=0
            widget_control, box_button_base, sensitive=0
            widget_control, i_boxbase, sensitive=0
            curbox=1
            totbox=1
            widget_control, curbox_field, set_value=strcompress(1)
            update_buttons
            calc_pixel, box_x(curbox), box_y(curbox)
            choose_standard_filter
            setup_ech_filter
            if mask_flag then erase_bad
            find_corners
            draw_echellogram
            draw_boxes
        end
        2:begin    ; Object Only
            widget_control, obj_name, sensitive=1
            widget_control, star_name, sensitive=0
            widget_control, itime, sensitive=1
            widget_control, stime, sensitive=0
            widget_control, coadds, sensitive=1
            widget_control, scoadds, sensitive=0
;            widget_control, boxnum, sensitive=1
            widget_control, box_button_base, sensitive=1
            widget_control, i_boxbase, sensitive=1
        end
        3:begin   ; Bright Object Only
            widget_control, obj_name, sensitive=1
            widget_control, star_name, sensitive=0
            widget_control, itime, sensitive=1
            widget_control, stime, sensitive=0
            widget_control, coadds, sensitive=1
            widget_control, scoadds, sensitive=0
;            widget_control, boxnum, sensitive=1
            widget_control, box_button_base, sensitive=1
            widget_control, i_boxbase, sensitive=1
        end
        4:begin      ; Not used
            widget_control, obj_name, sensitive=1
            widget_control, star_name, sensitive=1
            widget_control, itime, sensitive=1
            widget_control, stime, sensitive=1
            widget_control, coadds, sensitive=1
            widget_control, scoadds, sensitive=1
;            widget_control, boxnum, sensitive=0
            widget_control, box_button_base, sensitive=0
            widget_control, i_boxbase, sensitive=0
            curbox=1
            totbox=1
;            widget_control, boxnum, set_value=1
            widget_control, curbox_field, set_value=strcompress(1)
            update_buttons
            calc_pixel, box_x(curbox), box_y(curbox)
            choose_standard_filter
            setup_ech_filter
            if mask_flag then erase_bad
            find_corners
            draw_echellogram
            draw_boxes
        end
        5:begin        ; Not used
            widget_control, obj_name, sensitive=1
            widget_control, star_name, sensitive=1
            widget_control, itime, sensitive=1
            widget_control, stime, sensitive=1
            widget_control, coadds, sensitive=1
            widget_control, scoadds, sensitive=1
;            widget_control, boxnum, sensitive=0
            widget_control, box_button_base, sensitive=0
            widget_control, i_boxbase, sensitive=0
            curbox=1
            totbox=1
;            widget_control, boxnum, set_value=1
            widget_control, curbox_field, set_value=strcompress(1)
            update_buttons
            calc_pixel, box_x(curbox), box_y(curbox)
            choose_standard_filter
            setup_ech_filter
            if mask_flag then erase_bad
            find_corners
            draw_echellogram
            draw_boxes
        end
    endcase

endif

return
bad: on_ioerror, null
if ( unit gt 0 ) then free_lun, unit
err = widget_message("IO Error has occurred reading script.")
print, "IO Error has occured reading in script."
widget_control, set_value = "Error reading in script.", status_text
curbox = 1
totbox = 1
boxes(curbox) = {position, echtheta, gratheta, 'high res', 1, 0, $
                 filterarr(0).position, filterarr(0).wheel, $
                 'default', 10.0, 'default', 5.0, 92.0, 1, 1, fltarr(4), fltarr(4)}
slitlen = slitlist(boxes(curbox).slit).length
update_buttons
widget_control, curbox_field, set_value=strcompress(1)
calc_pixel, box_x(curbox), box_y(curbox)
if mask_flag then erase_bad
choose_standard_filter
setup_ech_filter
find_corners
draw_echellogram
draw_boxes

END
;**********************************************************

;**********************************************************
PRO clear_overlay
common shared_1
common shared_2
common shared_3
common shared_4
common shared_5

    ; set all flags to zero so that echellogram and boxes are redrawn,
    ; the lines are not.                            
    neon_stat=0
    argon_stat=0
    krypton_stat=0
    xenon_stat=0
    oh_stat=0
    user_stat=0
    draw_echellogram
    draw_boxes

END
;**********************************************************

;**********************************************************
PRO input_proff
common shared_1
common shared_2
common shared_3
common shared_4
common shared_5
common shared_6

END
;**********************************************************

;**********************************************************
PRO USER_LINES_SETUP
common shared_1
common shared_2
common shared_3
common shared_4
common shared_5
common shared_6


if ( user_stat eq 0 ) then begin ; turned on, then set flag
    user_stat=1
    user_lines
    return
endif
if (user_stat eq 1) then begin ; turned off, so set flag to zero
    user_stat=0
    draw_echellogram
    draw_boxes
    return
endif

END
;**********************************************************

;**********************************************************
PRO NEON_LINES_SETUP
common shared_1
common shared_2
common shared_3
common shared_4
common shared_5
common shared_6

if ( neon_stat eq 0 ) then begin ; turned on, then set flag and redraw
    neon_stat=1
    draw_echellogram
    draw_boxes
    return
endif
if (neon_stat eq 1) then begin ; turned off, then set flag to zero and redraw
    neon_stat=0
    draw_echellogram
    draw_boxes
    return
endif


END
;**********************************************************

;**********************************************************
PRO ARGON_LINES_SETUP
common shared_1
common shared_2
common shared_3
common shared_4
common shared_5
common shared_6

if ( argon_stat eq 0 ) then begin ; turned on, then set flag and redraw
    argon_stat=1
    draw_echellogram
    draw_boxes
    return
endif
if (argon_stat eq 1) then begin ; turned off, then set flag to zero and redraw
    argon_stat=0
    draw_echellogram
    draw_boxes
    return
endif


END
;**********************************************************

;**********************************************************
PRO KRYPTON_LINES_SETUP
common shared_1
common shared_2
common shared_3
common shared_4
common shared_5
common shared_6

if ( krypton_stat eq 0 ) then begin ; turned on, then set flag and redraw
    krypton_stat=1
    draw_echellogram
    draw_boxes
    return
endif
if (krypton_stat eq 1) then begin ; turned off, so set flag to zero and redraw
    krypton_stat=0
    draw_echellogram
    draw_boxes
    return
endif


END
;**********************************************************

;**********************************************************
PRO XENON_LINES_SETUP
common shared_1
common shared_2
common shared_3
common shared_4
common shared_5
common shared_6

if ( xenon_stat eq 0 ) then begin ; turned on, then set flag and redraw
    xenon_stat=1
    draw_echellogram
    draw_boxes
    return
endif
if (xenon_stat eq 1) then begin ; turned off, then set flag to zero and redraw
    xenon_stat=0
    draw_echellogram
    draw_boxes
    return
endif


END
;**********************************************************

;**********************************************************
PRO oh_LINES_SETUP
common shared_1
common shared_2
common shared_3
common shared_4
common shared_5
common shared_6

if ( oh_stat eq 0 ) then begin ; turned on, then set flag and redraw
    oh_stat=1
    draw_echellogram
    draw_boxes
    return
endif
if (oh_stat eq 1) then begin ; turned off, then set flag to zero and redraw
    oh_stat=0
    draw_echellogram
    draw_boxes
    return
endif


END
;**********************************************************

;**********************************************************
PRO DRAW_LINES
common shared_1
common shared_2
common shared_3
common shared_4
common shared_5
common shared_6

widget_control, draw, get_value=index
wset, index

techphi = fltarr(1)
techgra = fltarr(1)
blazelam = echsigma*cos(echgamma)*2.0*sin(echtheta)

if ( boxes(curbox).spmode eq 'high res' ) then begin
    if ( neon_stat eq 1 ) then begin
        read_lines, 'data/neon.lst', neon_lines, comments, neon_num, err
        if err eq 0 then begin
            for i=0, (neon_num-1) do begin
                if ((neon_lines(i) ge (l1*1000)) and (neon_lines(i) le $
                          (l2*1000))) then begin
                    echorder = fix(1000*blazelam/neon_lines(i) +0.5)
                    techphi(0) = asin(neon_lines(i)*float(echorder)/(2000.0* $
                         echsigma*cos(echgamma)))
                    techgra(0) = asin(neon_lines(i)*float(graorder)/(2000.0* $
                         grasigma*cos(grahaang)))
                    oplot, techphi, techgra, psym=4, color=neon_col
                    xyouts, techphi, (techgra+.001), neon_lines(i), $
                         alignment=0.7, orientation=chipang, color=neon_col
                endif
            endfor
        endif
    endif

    if ( argon_stat eq 1 ) then begin
        read_lines, 'data/argon.lst', argon_lines, comments, argon_num, err
        if err eq 0 then begin
        for i=0, (argon_num-1) do begin
            if ( (argon_lines(i) ge (l1*1000)) and (argon_lines(i) le $
                                                    (l2*1000)) ) then begin
                echorder = fix(1000*blazelam/argon_lines(i) +0.5)
                techphi(0) = asin(argon_lines(i)*float(echorder)/(2000.0* $
                                  echsigma*cos(echgamma)))
                techgra(0) = asin(argon_lines(i)*float(graorder)/(2000.0* $
                                  grasigma*cos(grahaang)))
                oplot, techphi, techgra, psym=4, color=argon_col
                xyouts, techphi, (techgra+.001), argon_lines(i), $
                  alignment=0.7, orientation=chipang, color=argon_col
            endif
        endfor
    endif
endif

    if ( krypton_stat eq 1 ) then begin
        read_lines, 'data/krypton.lst', krypton_lines, comments, $
          kryton_num, err
        if err eq 0 then begin
        for i=0, (krypton_num-1) do begin
            if ( (krypton_lines(i) ge (l1*1000)) and (krypton_lines(i) $
                      le (l2*1000))) then begin
                echorder = fix(1000*blazelam/krypton_lines(i) +0.5)
                techphi(0) = asin(krypton_lines(i)*float(echorder)/(2000.0* $
                                  echsigma*cos(echgamma)))
                techgra(0) = asin(krypton_lines(i)*float(graorder)/(2000.0* $
                                  grasigma*cos(grahaang)))
                oplot, techphi, techgra, psym=4, color=krypton_col
                xyouts, techphi, (techgra+.001), krypton_lines(i), $
                  alignment=0.7, orientation=chipang, color=krypton_col
            endif
        endfor
        endif
    endif
    if ( xenon_stat eq 1 ) then begin
        read_lines, 'data/xenon.lst', xenon_lines, comments, xenon_num, err
        if err eq 0 then begin
        for i=0, (xenon_num-1) do begin
            if ( (xenon_lines(i) ge (l1*1000)) and (xenon_lines(i) le $
                 (l2*1000)) ) then begin
                echorder = fix(1000*blazelam/xenon_lines(i) +0.5)
                techphi(0) = asin(xenon_lines(i)*float(echorder)/(2000.0* $
                                  echsigma*cos(echgamma)))
                techgra(0) = asin(xenon_lines(i)*float(graorder)/(2000.0* $
                                  grasigma*cos(grahaang)))
                oplot, techphi, techgra, psym=4, color=xenon_col
                xyouts, techphi, (techgra+.001), xenon_lines(i), $
                  alignment=0.7, orientation=chipang, color=xenon_col
            endif
        endfor
        endif
    endif
    if ( oh_stat eq 1 ) then begin
        read_lines, 'data/oh.lst', oh_lines, comments, oh_num, err
        if err eq 0 then begin
        for i=0, (oh_num-1) do begin
            if ( (oh_lines(i) ge (l1*1000)) and (oh_lines(i) le (l2*1000)) ) $
              then begin
                echorder = fix(1000*blazelam/oh_lines(i) +0.5)
                techphi(0) = asin(oh_lines(i)*float(echorder)/(2000.0* $
                                  echsigma*cos(echgamma)))
                techgra(0) = asin(oh_lines(i)*float(graorder)/(2000.0* $
                                  grasigma*cos(grahaang)))
                oplot, techphi, techgra, psym=4, color=oh_col
                xyouts, techphi, (techgra+.001), oh_lines(i), $
                  alignment=0.7, orientation=chipang, color=oh_col
            endif
        endfor
        endif
    endif

    if ( user_stat eq 1 ) then begin
            for i=0, (user_list_length-1) do begin
                if ((shifted_list(i) ge (l1*1000)) and (shifted_list(i) le $
                  (l2*1000))) then begin
                    echorder = fix(1000*blazelam/shifted_list(i) +0.5)
                    techphi(0) = asin(shifted_list(i)*float(echorder)/ $
                                      (2000.0*echsigma*cos(echgamma)))
                    techgra(0) = asin(shifted_list(i)*float(graorder)/ $
                                      (2000.0*grasigma*cos(grahaang)))
                    oplot, techphi, techgra, psym=4, color=user_col
                    xyouts, techphi, (techgra+.001), shifted_list(i), $
                         alignment=0.7, orientation=chipang, color=user_col
                endif
            endfor
    endif

endif

if ( boxes(curbox).spmode eq 'low res' ) then begin
    if ( neon_stat eq 1 ) then begin
        read_lines, 'data/neon.lst', neon_lines, comments, neon_num, err
        if err eq 0 then begin
        for i=0, (neon_num-1) do begin
            if ( (neon_lines(i) ge (l1*1000)) and (neon_lines(i) le $
                                                   (l2*1000)) ) then begin
                echorder = fix(1000*blazelam/neon_lines(i) +0.5)
                techphi(0) = asin(neon_lines(i)*float(echorder)/(2000.0* $
                                  echsigma*cos(echgamma)))
                techgra(0) = asin(neon_lines(i)*float(graorder)/(2000.0* $
                                  grasigma*cos(grahaang)))
                oplot, techphi, techgra, psym=4, color=neon_col
                xyouts, techphi, (techgra-.0002), neon_lines(i), $
                  alignment=0.28, color=neon_col
            endif
        endfor
        endif
    endif
    if ( argon_stat eq 1 ) then begin
        read_lines, 'data/argon.lst', argon_lines, comments, argon_num, err
        if err eq 0 then begin
        for i=0, (argon_num-1) do begin
            if ( (argon_lines(i) ge (l1*1000)) and (argon_lines(i) le $
               (l2*1000)) ) then begin
                echorder = fix(1000*blazelam/argon_lines(i) +0.5)
                techphi(0) = asin(argon_lines(i)*float(echorder)/(2000.0* $
                                  echsigma*cos(echgamma)))
                techgra(0) = asin(argon_lines(i)*float(graorder)/(2000.0* $
                                  grasigma*cos(grahaang)))
                oplot, techphi, techgra, psym=4, color=argon_col
                xyouts, techphi, (techgra-.0002), argon_lines(i), $
                  alignment=0.28, color=argon_col
            endif
        endfor
        endif
    endif
    if ( krypton_stat eq 1 ) then begin
        read_lines, 'data/krypton.lst', krypton_lines, comments, krypton_num, $
          err
        if err eq 0 then begin
        for i=0, (krypton_num-1) do begin
            if ( (krypton_lines(i) ge (l1*1000)) and (krypton_lines(i) le $
               (l2*1000)) ) then begin
                echorder = fix(1000*blazelam/krypton_lines(i) +0.5)
                techphi(0) = asin(krypton_lines(i)*float(echorder)/(2000.0* $
                                  echsigma*cos(echgamma)))
                techgra(0) = asin(krypton_lines(i)*float(graorder)/(2000.0* $
                                  grasigma*cos(grahaang)))
                oplot, techphi, techgra, psym=4, color=krypton_col
                xyouts, techphi, (techgra-.0002), krypton_lines(i), $
                  alignment=0.28, color=krypton_col
            endif
        endfor
        endif
    endif
    if ( xenon_stat eq 1 ) then begin
        read_lines, 'data/xenon.lst', xenon_lines, comments, xenon_num, err
        if err eq 0 then begin
        for i=0, (xenon_num-1) do begin
            if ( (xenon_lines(i) ge (l1*1000)) and (xenon_lines(i) le $
               (l2*1000)) ) then begin
                echorder = fix(1000*blazelam/xenon_lines(i) +0.5)
                techphi(0) = asin(xenon_lines(i)*float(echorder)/(2000.0* $
                                  echsigma*cos(echgamma)))
                techgra(0) = asin(xenon_lines(i)*float(graorder)/(2000.0* $
                                  grasigma*cos(grahaang)))
                oplot, techphi, techgra, psym=4, color=xenon_col
                xyouts, techphi, (techgra-.0002), xenon_lines(i), $
                  alignment=0.28, color=xenon_col
            endif
        endfor
        endif
    endif
    if ( oh_stat eq 1 ) then begin
        read_lines, 'data/oh.lst', oh_lines, comments, oh_num, err
        if err eq 0 then begin
        for i=0, (oh_num-1) do begin
            if ( (oh_lines(i) ge (l1*1000)) and (oh_lines(i) le $
              (l2*1000)) ) then begin
                echorder = fix(1000*blazelam/oh_lines(i) +0.5)
                techphi(0) = asin(oh_lines(i)*float(echorder)/(2000.0* $
                                  echsigma*cos(echgamma)))
                techgra(0) = asin(oh_lines(i)*float(graorder)/(2000.0* $
                                  grasigma*cos(grahaang)))
                oplot, techphi, techgra, psym=4, color=oh_col
                xyouts, techphi, (techgra-.0002), oh_lines(i), $
                  alignment=0.28, color=oh_col
            endif
        endfor
        endif
    endif
    
    if ( user_stat eq 1 ) then begin
        for i=0, (user_list_length-1) do begin
            if ( (shifted_list(i) ge (l1*1000)) and (shifted_list(i) le $
              (l2*1000)) ) then begin
                echorder = fix(1000*blazelam/shifted_list(i) +0.5)
                techphi(0) = asin(shifted_list(i)*float(echorder)/(2000.0* $
                                  echsigma*cos(echgamma)))
                techgra(0) = asin(shifted_list(i)*float(graorder)/(2000.0* $
                                  grasigma*cos(grahaang)))
                oplot, techphi, techgra, psym=4, color=user_col
                xyouts, techphi, (techgra-.0002), shifted_list(i), $
                  alignment=0.28, color=user_col
            endif
        endfor
    endif
endif

END
;**********************************************************


;**********************************************************
PRO INPUT_WAVELENGTH_EVENT, Event
common shared_1
common shared_2
common shared_3
common shared_4
common shared_5
common shared_6
common shared_ao

END
;**********************************************************

;**********************************************************
PRO INPUT_WAVELENGTH
common shared_1
common shared_2
common shared_3
common shared_4
common shared_5
common shared_6
common shared_ao

END
;**********************************************************

;**********************************************************
PRO SETUP_ECH_FILTER
common shared_1
common shared_2
common shared_3
common shared_4
common shared_5
common shared_6
common shared_ao

if filtermenu_index(curbox) ne 18 then begin

if ( boxes(curbox).spmode eq 'high res' ) then begin
    blazelam = echsigma*cos(echgamma)*2.0*sin(echtheta)
;	leftlam:	array of leftmost wavelengths for plotted orders
;	rightlam:	array of rightmost wavelengths for plotted orders.
;	leftech:	array of echelle angles for leftmost wavelengths.
;	rightech:	array of echelle angles for rightmost wavelengths.
;	leftgra:	array of grating angles for leftmost wavelengths.
;	rightgra:	array of grating angles for rightmost wavelengths.
;	numorders:	number of orders for current filter
    order2 = fix((blazelam/l1) + 0.5)
    order1 = fix((blazelam/l2) + 0.5)
    numorders = order2 - order1
    slitlen = slitlist(boxes(curbox).slit).length
    
    chiprot = !pi*chipang/180.0
    st = sin(chiprot)
    ct = cos(chiprot)
; calculate the middle angles for the crossover wavelengths
    for i = order1, order2 do begin
        echorder = i
        leftlam(i) = blazelam/(float(echorder)+0.5)
        rightlam(i) = blazelam/(float(echorder)-0.5)
        if ( leftlam(i) lt l1 ) then leftlam(i) = l1
        if ( rightlam(i) gt l2 ) then rightlam(i) = l2
        leftech(i) = asin(leftlam(i)*float(echorder)/(2.0*echsigma* $
                                                      cos(echgamma)))
        leftgra(i) = asin(leftlam(i)*float(graorder)/(2.0*grasigma* $
                                                      cos(grahaang)))
        rightech(i) = asin(rightlam(i)*float(echorder)/(2.0*echsigma* $
                                                        cos(echgamma)))
        rightgra(i) = asin(rightlam(i)*float(graorder)/(2.0*grasigma* $
                                                        cos(grahaang)))
        uleftx(i) = 0.
        ulefty(i) = slitlen/yarc/2.
        urightx(i) = 0.
        urighty(i) = slitlen/yarc/2.
        tx = uleftx(i) * ct - ulefty(i)*st
        ty = uleftx(i)*st + ulefty(i)*ct
        uleftx(i) = tx
        ulefty(i) = ty
        tx = urightx(i) * ct - urighty(i)*st
        ty = urightx(i)*st + urighty(i)*ct
        urightx(i) = tx
        urighty(i) = ty
        echdel = atan( uleftx(i) * pix / few )
        gradel = atan( ulefty(i) * pix / feh )
        
        uleftx(i) = asin( (sin(leftech(i)) + $
                           sin( (leftech(i) + echdel) )) / 2.0 )
        ulefty(i) = asin( (sin(leftgra(i)+grahaang) + $
                           sin( (leftgra(i) -grahaang + gradel) )) / (2.0* $
                           cos(grahaang)) )
        echdel = atan( urightx(i) * pix / few )
        gradel = atan( urighty(i) * pix / feh )
        urightx(i) = asin( (sin(rightech(i)) + $
                            sin( (rightech(i) + echdel) )) / 2.0 )
        urighty(i) = asin( (sin(rightgra(i)+grahaang) + $
                            sin( (rightgra(i) -grahaang + gradel) )) / (2.0* $
                            cos(grahaang)) )
        
        
        
        lleftx(i) = 0.
        llefty(i) = 0. - (slitlen/yarc/2.)
        lrightx(i) = 0.
        lrighty(i) = 0. - (slitlen/yarc/2.)
        tx = lleftx(i) * ct - llefty(i)*st
        ty = lleftx(i)*st + llefty(i)*ct
        lleftx(i) = tx
        llefty(i) = ty
        tx = lrightx(i) * ct - lrighty(i)*st
        ty = lrightx(i)*st + lrighty(i)*ct
        lrightx(i) = tx
        lrighty(i) = ty
        echdel = atan( lleftx(i) * pix / few )
        gradel = atan( llefty(i) * pix / feh )
        lleftx(i) = asin( (sin(leftech(i)) + $
                           sin( (leftech(i) + echdel) )) / 2.0 )
        llefty(i) = asin( (sin(leftgra(i)+grahaang) + $
                           sin( (leftgra(i) -grahaang + gradel) )) / (2.0* $
                                                                      cos(grahaang)) )
        echdel = atan( lrightx(i) * pix / few )
        gradel = atan( lrighty(i) * pix / feh )
        lrightx(i) = asin( (sin(rightech(i)) + $
                            sin( (rightech(i) + echdel) )) / 2.0 )
        lrighty(i) = asin( (sin(rightgra(i)+grahaang) + $
                            sin( (rightgra(i) -grahaang + gradel) )) / (2.0* $
                                                                        cos(grahaang)) )
    endfor
        
    
; Calculate positions of 1000 wavelengths within band.
    dellam = (l2-l1)/999.0
    for i = 0, 999 do begin
        dlambda(i) = l1 + dellam*float(i)
        tmp = (blazelam/dlambda(i)) + 0.5
        echorder = fix(tmp)
        echphi(i) = asin(dlambda(i)*float(echorder)/(2.0*echsigma* $
                                                     cos(echgamma)))
        graphi(i) = asin(dlambda(i)*float(graorder)/(2.0*grasigma* $
                                                     cos(grahaang)))
        lgraphi(i) = asin(dlambda(i)*float(graorder-1)/(2.0*grasigma* $
                                                        cos(grahaang)))
        hgraphi(i) = asin(dlambda(i)*float(graorder+1)/(2.0*grasigma* $
                                                        cos(grahaang)))
    endfor
; Calculate plotting range.
; Start with width of echellogram for about 1.0 microns.
    data_x1 = 1.08529
    data_x2 = 1.11490
    data_y1 = graphi(0)
    data_y2 = graphi(0)
    for i = 0, 999 do begin
        if ( echphi(i) le data_x1) then data_x1 = echphi(i)
        if ( echphi(i) ge data_x2) then data_x2 = echphi(i)
        if ( graphi(i) le data_y1) then data_y1 = graphi(i)
        if ( graphi(i) ge data_y2) then data_y2 = graphi(i)
;	if ( lgraphi(i) le data_y1) then data_y1 = lgraphi(i)
;	if ( lgraphi(i) ge data_y2) then data_y2 = lgraphi(i)
;	if ( hgraphi(i) le data_y1) then data_y1 = hgraphi(i)
;	if ( hgraphi(i) ge data_y2) then data_y2 = hgraphi(i)
    endfor
;    xrange = (data_x2 - data_x1)*0.7
    xrange = (data_x2 - data_x1)*0.75
    xmid = (data_x2 + data_x1)/2.0
;    yrange = (data_y2 - data_y1)*0.7
    yrange = (data_y2 - data_y1)*0.75
    ymid = (data_y2 + data_y1)/2.0
    if ( aspect gt 1.0 ) then yrange = yrange * aspect
    if ( aspect le 1.0 ) then xrange = xrange / aspect
    data_x1 = xmid - xrange
    data_x2 = xmid + xrange
    data_y1 = ymid - yrange
    data_y2 = ymid + yrange
endif
if ( boxes(curbox).spmode eq 'low res' ) then begin
;	leftgra:	bottommost wavelengths.
;	rightgra:	topmost wavelengths.
;  order2 = fix((blazelam/l1) + 0.5)
;  order1 = fix((blazelam/l2) + 0.5)
;  numorders = order2 - order1
    slitlen = slitlist(boxes(curbox).slit).width
    chiprot = !pi*chipang/180.0
    st = sin(chiprot)
    ct = cos(chiprot)
; calculate the middle angles for the crossover wavelengths
    leftgra(1) = asin(l1*float(graorder)/(2.0*grasigma* $
                                          cos(grahaang)))
    rightgra(1) = asin(l2*float(graorder)/(2.0*grasigma* $
                                           cos(grahaang)))
    uleftx(1) = slitlen/yarc/2.
    urightx(1) = 3.0*slitlen/yarc/2.
    lleftx(1) = slitlen/yarc/2.
    lrightx(1) = 3.0*slitlen/yarc/2.
;    uleftx(1) = 0. - slitlen/yarc/2.
;    urightx(1) = slitlen/yarc/2.
;    lleftx(1) = 0. - slitlen/yarc/2.
;    lrightx(1) = slitlen/yarc/2.
    ulefty(1) = rightgra(1)
    urighty(1) = rightgra(1)
    llefty(1) = leftgra(1)
    lrighty(1) = leftgra(1)
    
; Calculate positions of 1000 wavelengths within band.
    dellam = (l2-l1)/999.0
    for i = 0, 999 do begin
        dlambda(i) = l1 + dellam*float(i)
        graphi(i) = asin(dlambda(i)*float(graorder)/(2.0*grasigma* $
                                                     cos(grahaang)))
    endfor
    
    yrange = (urighty(1) - lrighty(1))/1.5
    xrange = 700.
    if ( yaspect gt 1.0 ) then yrange = yrange * yaspect
    if ( yaspect lt 1.0 ) then xrange = xrange / yaspect
    ymid = (urighty(1) + lrighty(1))/2.0
    data_x1 = 0.0 - xrange
    data_x2 = 0.0 + xrange
    data_y1 = ymid - yrange
    data_y2 = ymid + yrange
    yrange = (urighty(1) - lrighty(1))/1.8
    ymid = (urighty(1) + lrighty(1))/2.0
    if ( (yaspect gt 3) and (yaspect le 8) ) then begin
    	llefty(2) = ymid - (yrange*yaspect/3.5)
    	ulefty(2) = ymid + yrange
    endif
    if ( yaspect gt 8 ) then begin
    	llefty(2) = ymid - (yrange*yaspect/7)
    	ulefty(2) = ymid + yrange*1.5
    endif
    if ( yaspect le 3 ) then begin
    	llefty(2) = ymid - yrange
	ulefty(2) = ymid + (yrange*.95)
    endif
endif

endif

END
;**********************************************************

;**********************************************************
PRO find_corners

common shared_1
common shared_2
common shared_3
common shared_4
common shared_5
common shared_6
common shared_ao

; Determine corners of current box in grating angles.
; First rotate the box corners 4.657 degrees

; MODIFIED  GWD  Nov. 2018 to replace 1024 with 2048 pixel array

;boxes(curbox).x(0) = -512.0
;boxes(curbox).y(0) = -512.0
;boxes(curbox).x(1) = 512.0
;boxes(curbox).y(1) = -512.0
;boxes(curbox).x(2) = 512.0
;boxes(curbox).y(2) = 512.0
;boxes(curbox).x(3) = -512.0
;boxes(curbox).y(3) = 512.0
boxes(curbox).x(0) = -1024.0
boxes(curbox).y(0) = -1024.0
boxes(curbox).x(1) = 1024.0
boxes(curbox).y(1) = -1024.0
boxes(curbox).x(2) = 1024.0
boxes(curbox).y(2) = 1024.0
boxes(curbox).x(3) = -1024.0
boxes(curbox).y(3) = 1024.0

chiprot = !pi*chipang/180.0

st = sin(chiprot)
ct = cos(chiprot)
for i= 0, 3 do begin
    tx = boxes(curbox).x(i) *ct - boxes(curbox).y(i)*st
    ty = boxes(curbox).x(i) *st + boxes(curbox).y(i)*ct
    boxes(curbox).x(i) = tx
    boxes(curbox).y(i) = ty
endfor

if ( boxes(curbox).spmode eq 'high res' ) then begin
; Now final delta betas for given pixel locations and calculate
; equivalent phi's for plotting
    for i = 0, 3 do begin
        echdel = atan( boxes(curbox).x(i) * pix / dist_few(i) )
        gradel = atan( boxes(curbox).y(i) * pix / dist_feh(i) )
        boxes(curbox).x(i) = asin( (sin(boxes(curbox).echphi) + $
                                   sin( (boxes(curbox).echphi + echdel) $
                                       )) / 2.0 )
        boxes(curbox).y(i) = asin( (sin(boxes(curbox).graphi+grahaang) + $
                                    sin( (boxes(curbox).graphi -grahaang + $
                                          gradel) )) / (2.0* $
                                                        cos(grahaang)) )
    endfor
endif
if ( boxes(curbox).spmode eq 'low res' ) then begin
    for i = 0, 3 do begin
        gradel = atan( boxes(curbox).y(i) * pix / dist_feh(i) )
        boxes(curbox).y(i) = asin( (sin(boxes(curbox).graphi+grahaang) + $
                                    sin( (boxes(curbox).graphi -grahaang + $
                                          gradel) )) / (2.0* $
                                                        cos(grahaang)) )
    endfor
endif

END
;**********************************************************

;**********************************************************
PRO CHOOSE_STANDARD_FILTER

common shared_1
common shared_2
common shared_3
common shared_4
common shared_5
common shared_6
common shared_ao

;for i = 0, (numfil-1) do begin
;    if ( (boxes(curbox).filter eq filterarr(i).position) and $
;       (boxes(curbox).wheel eq filterarr(i).wheel) ) then begin
i=filtermenu_index(curbox)
        l1 = filterarr(i).l1
        l2 = filterarr(i).l2
        graorder = filterarr(i).graorder
        aspect = filterarr(i).aspect
        yaspect = filterarr(i).yaspect
;    endif
;endfor

END
;**********************************************************

;**********************************************************
; this routine was added by JLW on 2/10/99
; it is called by draw_event to update the x and y pixel
; locations and the lambda and ech order as well.
; this routine was internal to draw_event, but i extracted
; it so that it is called when arrow keys are hit to move
; the box.

pro calc_pixel, event_x, event_y
common shared_1
common shared_2
common shared_3
common shared_4
common shared_5
common shared_6


if filtermenu_index(curbox) ne 18 then begin

blazelam = echsigma*cos(echgamma)*2.0*sin(echtheta)
if boxes(curbox).spmode eq 'high res' then begin
        echelle = (data_x2-data_x1)*event_x/scrx + data_x1
        grating = (data_y2-data_y1)*event_y/scry + data_y1
        lam = 2.0*grasigma*cos(grahaang)*sin(grating)/float(graorder)
;  MODIFIED GWD Nov 2018
;       ypix = (((1000.0*lam * graorder)-(24.07*sin(boxes(curbox).graphi))) $
;		/  0.00084374) + 512.5
        ypix = (((1000.0*lam * graorder)-(24.07*sin(boxes(curbox).graphi))) $
		/  0.0005625) + 1024.5
        lechorder = fix(blazelam/lam)
        uechorder = lechorder + 1
        llam = sin(echelle)*2.0*echsigma*cos(echgamma)/float(lechorder)
        ulam = sin(echelle)*2.0*echsigma*cos(echgamma)/float(uechorder)
        ldiff = abs(llam - lam)
        udiff = abs(ulam - lam)
        lam = llam*1000.0
        ord = lechorder
        if ( ldiff gt udiff ) then begin
            lam = ulam*1000.0
            ord = uechorder
        endif
        widget_control, echbox, $
		set_value=string(round2(boxes(curbox).echphi*180.0/!pi))
        widget_control, grabox, set_value=$
		string(round2((boxes(curbox).graphi+ grahaang)*180.0/!pi))
        widget_control, lambox, set_value=string(lam)
        widget_control, ordbox, set_value=string(ord)
 
; MODIFIED GWD Nov 2018
;        xpix = (((lam*ord)-(85.5034*sin(boxes(curbox).echphi)))/0.00110235) $
;          + 512.5
        xpix = (((lam*ord)-(85.5034*sin(boxes(curbox).echphi)))/0.0007349) $
          + 1024.5

        widget_control, xpixbox, set_value=string(xpix)
        widget_control, ypixbox, set_value=string(ypix)
endif
if boxes(curbox).spmode eq 'low res' then begin
        grating = (data_y2-data_y1)*event_y/scry + data_y1
        lam = 1000.0*2.0*grasigma*cos(grahaang)*sin(grating)/float(graorder)

; MODIFIED GWD Nov 2018        
;        ypix = (((lam * graorder)-(24.07*sin(boxes(curbox).graphi)))/  $
;                0.00084374) + 512.5
;        xpix = 512
        ypix = (((lam * graorder)-(24.07*sin(boxes(curbox).graphi)))/  $
                0.0005625) + 1024.5
        xpix = 1024
        widget_control, echbox, set_value='179.7'
        widget_control, grabox, set_value=$
		string(round2((boxes(curbox).graphi + grahaang)*180.0/!pi))
        widget_control, lambox, set_value=string(lam)
        widget_control, ordbox, set_value='None'
        widget_control, xpixbox, set_value=string(xpix)
        widget_control, ypixbox, set_value=string(ypix)
endif
if boxes(curbox).spmode eq 'imaging' then begin
        widget_control, echbox, set_value='179.7'
        widget_control, grabox, set_value='0th order'
        widget_control, lambox, set_value='Undispersed'
        widget_control, ordbox, set_value='None'
;        xpix = 512
;        ypix = 512
        xpix = 1024
        ypix = 1024
        widget_control, xpixbox, set_value=string(xpix)
        widget_control, ypixbox, set_value=string(ypix)
endif

endif

end
;********************************************************************

;*********************************************************************
PRO DRAW_EVENT, event
common shared_1
common shared_2
common shared_3
common shared_4
common shared_5
common shared_6
common shared_ao

widget_control, draw, get_value=index
wset, index

if tag_names(event, /structure_name) eq 'WIDGET_TRACKING' then begin
    in_draw=event.enter
endif else begin

csr_x=event.x
csr_y=event.y

if blank then begin
    if event.type eq 0 then begin
        valid=0
        not_in_box=1
        for i=1, totbox do begin
            if (check_in_box(event.x, event.y, [box_x(i)-10, box_x(i)+10, $
              box_x(i)+10, box_x(i)-10], [box_y(i)-15, box_y(i)-15, $
              box_y(i)+15, box_y(i)+15], 'imaging') ne 0) and $
              (filtermenu_index(i) eq current_fil) and $
              (boxes(i).spmode eq sptype) then begin
                valid=i
                not_in_box=0
                if i eq curbox then i=totbox+1
            endif
        endfor

        if not_in_box then begin
            box_x(curbox)=event.x
            box_y(curbox)=event.y
            ; fill boxes(curbox) with new values
            boxes(curbox).spmode=sptype
            if filtermenu_index(curbox) ne current_fil then $
              boxes(curbox).blocker=2
            filtermenu_index(curbox)=current_fil
; Changed temporarily for first light to 179.99 degrees
    boxes(curbox).echphi = roundphi(179.7*!dtor)
;            boxes(curbox).echphi = !pi
            boxes(curbox).graphi = 0.0
            ; Set the boxes filter and wheel parameters
            boxes(curbox).filter = filterarr(current_fil).position
            boxes(curbox).wheel = filterarr(current_fil).wheel
            boxes(curbox).slit=current_slit
            slitlen = slitlist(boxes(curbox).slit).length
        endif

        if valid gt 0 then curbox=valid
        update_buttons
        draw_echellogram
        draw_boxes

    endif
endif else begin  ; not blank

if sptype ne 'imaging' then begin
        calc_pixel, event.x, event.y


 if ( event.type eq 0 ) then begin ; Button pressed
	numvalid=0
	not_in_box=1
	valid_boxes=intarr(maxbox)
        ; set window since cw_toggle messes up data coords
        widget_control, draw, get_value=index
        wset, index
	dat_xy=convert_coord(event.x, event.y, /device, /to_data)
	for i=1, totbox do begin
		if (check_in_box(event.x, event.y, boxes(i).x, $
			boxes(i).y, sptype) ne 0) and (filtermenu_index(i) $ 
			eq current_fil) and (boxes(i).spmode eq sptype) $
			then begin
			valid_boxes(numvalid)=i
			numvalid=numvalid+1
			not_in_box=0
			if (i) eq curbox then begin
				xoff=event.x-box_x(i)
				yoff=event.y-box_y(i)
        			moving = 1
				i=totbox
			endif			
		endif
	endfor

	if (not_in_box) then begin
		; fill boxes(curbox) with new values
		if filtermenu_index(curbox) ne current_fil then $
                  boxes(curbox).blocker=2
                filtermenu_index(curbox)=current_fil
		boxes(curbox).spmode = sptype
		if boxes(curbox).spmode eq 'high res' THEN BEGIN
			echelle = (data_x2-data_x1)*event.x/scrx + data_x1
        		boxes(curbox).echphi = roundphi(echelle)
	    	endif
     		grating = (data_y2-data_y1)*event.y/scry + data_y1
        	boxes(curbox).graphi = roundphi(grating)
		; Set the boxes filter and wheel parameters
		boxes(curbox).filter = filterarr(current_fil).position
		boxes(curbox).wheel = filterarr(current_fil).wheel
		boxes(curbox).slit=current_slit
		slitlen = slitlist(boxes(curbox).slit).length
		calc_pixel, box_x(curbox), box_y(curbox)
		update_buttons
		find_corners
;		erase_curbox
		draw_echellogram
		draw_boxes
		xoff=0
		yoff=0
    		moving = 1
	endif

	if (numvalid gt 0 and moving eq 0) then begin
		curbox=max(valid_boxes)
		xoff=event.x-box_x(curbox)
		yoff=event.y-box_y(curbox)
    		moving = 1
	    slitlen = slitlist(boxes(curbox).slit).length
	    update_buttons
	    calc_pixel, event.x, event.y
 	    choose_standard_filter
 	    setup_ech_filter
	    if mask_flag then erase_bad
	    find_corners
	    draw_echellogram
	    draw_boxes

	endif
 ENDIF

 if ( event.type eq 1 ) and (moving eq 1)then begin ; Button released
	box_x(curbox)=event.x-xoff
	box_y(curbox)=event.y-yoff
	calc_pixel, event.x, event.y
        moving = 0
 ENDIF

 if ( (event.type eq 2) and (moving eq 1) ) then begin ; moving with button

	mov_x=event.x-xoff
	mov_y=event.y-yoff
    if boxes(curbox).spmode eq 'high res' THEN BEGIN
        echelle = (data_x2-data_x1)*mov_x/scrx + data_x1
        boxes(curbox).echphi = roundphi(echelle)
    endif
        grating = (data_y2-data_y1)*mov_y/scry + data_y1
        boxes(curbox).graphi = roundphi(grating)
;	calc_pixel, mov_x, mov_y
;	calc_pixel, event.x, event.y
        erase_curbox            ; Erase current box
        find_corners
        erase_curbox            ; Draw current box
        if mask_flag then begin
            erase_bad
            find_bad
            erase_bad
        endif
 ENDIF

endif else begin  ; in imaging
    if event.type eq 0 then begin
        valid=0
        not_in_box=1
        for i=1, totbox do begin
            if (check_in_box(event.x, event.y, [box_x(i)-10, box_x(i)+10, $
              box_x(i)+10, box_x(i)-10], [box_y(i)-15, box_y(i)-15, $
              box_y(i)+15, box_y(i)+15], sptype) ne 0) and $
              (filtermenu_index(i) eq current_fil) and $
              (boxes(i).spmode eq sptype) then begin
                valid=i
                not_in_box=0
                if i eq curbox then i=totbox+1
            endif
        endfor

        if not_in_box then begin
            widget_control, draw, get_value=index
            wset, index 
            if curbox eq 0 then curbox=old_curbox
            box_x(curbox)=event.x
            box_y(curbox)=event.y
            ; fill boxes(curbox) with new values
            boxes(curbox).spmode='imaging'
            if filtermenu_index(curbox) ne current_fil then $
              boxes(curbox).blocker=2
            filtermenu_index(curbox)=current_fil
; Changed temporarily for first light to 179.99 degrees
    boxes(curbox).echphi = roundphi(179.7*!dtor)
;            boxes(curbox).echphi = !pi
            boxes(curbox).graphi = 0.0
            ; Set the boxes filter and wheel parameters
            boxes(curbox).filter = filterarr(current_fil).position
            boxes(curbox).wheel = filterarr(current_fil).wheel
            boxes(curbox).slit=current_slit
            slitlen = slitlist(boxes(curbox).slit).length
        endif

        if valid gt 0 then curbox=valid
        update_buttons
        draw_echellogram
        draw_boxes

    endif
endelse

endelse

endelse

END
;**********************************************************

;**********************************************************
PRO draw_echellogram

common shared_1
common shared_2
common shared_3
common shared_4
common shared_5
common shared_6
common shared_ao

widget_control, draw, get_value=index
wset, index

if (filterarr(filtermenu_index(curbox)).name eq 'BLANK') then BEGIN
    plot, [1,46], [1,46], xtitle='arcseconds', $
      ytitle='arcseconds', xstyle=1, ystyle=1, /nodata
    xyouts, 23, 23, 'BLANK', /data, charsize=8, align=0.5, $
      charthick=10
Endif else begin

if (boxes(curbox).spmode eq 'high res') THEN BEGIN
    plot, echphi, graphi, xrange=[data_x1, data_x2],$
      yrange=[data_y1, data_y2], color=0, $
      position=[0,0, 1, 1], xstyle=1, ystyle=1
    templam = 1000.0 * leftlam
    xyouts, lleftx(order1:order2), llefty(order1:order2), $
      templam(order1:order2), alignment=1.1, orientation=chipang, $
      charsize=1.2, color=text_col
    templam = 1000.0 * rightlam
    xyouts, lrightx(order1:order2), lrighty(order1:order2), $
      templam(order1:order2), alignment=0.3, orientation=chipang, $
      charsize=1.2, color=text_col
    for i = order1, order2 do begin
        pxval = [lleftx(i),lrightx(i),urightx(i),uleftx(i), lleftx(i)]
        pyval = [llefty(i),lrighty(i),urighty(i),ulefty(i), llefty(i)]
        polyfill, pxval, pyval, color=spec_col
        oplot, pxval, pyval, color=border_col
    endfor
;    oplot, echphi, lgraphi, psym=4, color=border_col
;    oplot, echphi, hgraphi, psym=4, color=border_col
    ; since echelle is redrawn, call fill_curbox_arrays to update box_x
    ; and box_y with new data/device coordinate conversion
    fill_curbox_arrays

endif
if (boxes(curbox).spmode eq 'low res') then BEGIN
    plot, echphi, graphi, xrange=[data_x1, data_x2],$
      yrange=[data_y1, data_y2], color=0, $
      position=[0,0, 1, 1], xstyle=1, ystyle=1
    pxval = [lleftx(1),lrightx(1),urightx(1),uleftx(1), lleftx(1)]
    pyval = [llefty(1),lrighty(1),urighty(1),ulefty(1), llefty(1)]
    templam = l1 * 1000.0
    xyouts, 0, llefty(2), templam, charsize=1.2, alignment=0.7
    templam = l2 * 1000.0
    xyouts, 0, ulefty(2), templam, charsize=1.2, alignment=0.7
    polyfill, pxval, pyval, color=spec_col
    oplot, pxval, pyval, color=border_col
    ; since echelle is redrawn, call fill_curbox_arrays to update box_x
    ; and box_y with new data/device coordinate conversion
    fill_curbox_arrays

ENDif
if (boxes(curbox).spmode eq 'imaging') then BEGIN
    plot, [1,46], [1,46], xtitle='arcseconds', $
      ytitle='arcseconds', xstyle=1, ystyle=1, /nodata
    xyouts, 23, 23, 'IMAGING', /data, charsize=8, align=0.5, $
      charthick=10
    boxes(curbox).spmode = 'imaging'

ENDif
draw_lines

endelse


END
;**********************************************************

;**********************************************************
PRO draw_boxes
common shared_1
common shared_2
common shared_3
common shared_4
common shared_5
common shared_6
common shared_ao

widget_control, draw, get_value=index
wset, index

for i = 1, totbox do begin
    if ( boxes(i).spmode eq boxes(curbox).spmode ) then begin
	if filtermenu_index(i) eq current_fil then begin 
            if filtermenu_index(i) ne 18 then begin
             if ( boxes(i).spmode eq 'high res' ) then begin
                device, set_graphics = 6
                thick = 3
                col = off_col
                if ( i eq curbox ) then col = sel_col
                plots, [boxes(i).x(0),boxes(i).x(1)], $
                  [boxes(i).y(0),boxes(i).y(1)], $
                  color = col, thick = thick
                plots, [boxes(i).x(2),boxes(i).x(1)], $
                  [boxes(i).y(2),boxes(i).y(1)], $
                  color = col, thick = thick
                plots, [boxes(i).x(2),boxes(i).x(3)], $
                  [boxes(i).y(2),boxes(i).y(3)], $
                  color = col, thick = thick
                plots, [boxes(i).x(0),boxes(i).x(3)], $
                  [boxes(i).y(0),boxes(i).y(3)], $
                  color = col, thick = thick
                boxlab = strtrim(string(i),1)
                xlab = (0.9 * boxes(i).x(0)) + 0.1*boxes(i).x(1)
                ylab = (0.9 * boxes(i).y(0)) + 0.1*boxes(i).y(3)
                xyouts, xlab, ylab, boxlab, $
                  charsize = 4, color = col, charthick = 3
                device, set_graphics = 3
                if mask_flag then begin
                    if i eq curbox then begin
                        find_bad
                        erase_bad
                    endif
                endif
            endif
            if (boxes(curbox).spmode eq 'low res') then BEGIN
                device, set_graphics = 6
                thick = 3
                col = off_col
                if ( i eq curbox ) then col = sel_col
                plots, [boxes(i).x(0),boxes(i).x(1)], $
                  [boxes(i).y(0),boxes(i).y(1)], $
                  color = col, thick = thick
                plots, [boxes(i).x(2),boxes(i).x(1)], $
                  [boxes(i).y(2),boxes(i).y(1)], $
                  color = col, thick = thick
                plots, [boxes(i).x(2),boxes(i).x(3)], $
                  [boxes(i).y(2),boxes(i).y(3)], $
                  color = col, thick = thick
                plots, [boxes(i).x(0),boxes(i).x(3)], $
                  [boxes(i).y(0),boxes(i).y(3)], $
                  color = col, thick = thick
                boxlab = strtrim(string(i),1)
                xlab = (0.9 * boxes(i).x(0)) + 0.1*boxes(i).x(1)
                ylab = (0.9 * boxes(i).y(0)) + 0.1*boxes(i).y(3)
                xyouts, xlab, ylab, boxlab, $
                  charsize = 4, color = col, charthick = 3
                device, set_graphics = 3
                if mask_flag then begin
                    if i eq curbox then begin
                        find_bad
                        erase_bad
                    endif
                endif
            ENDif
            if (boxes(curbox).spmode eq 'imaging') then BEGIN
                device, set_graphics = 6
                thick = 3
                col = off_col
                if ( i eq curbox ) then col = sel_col
                xyouts, box_x(i), box_y(i)-15, strcompress(i, $
                  /remove_all), charthick=thick, color=col, charsize=4, $
                  /device, alignment=0.5
                device, set_graphics = 3

            ENDif
        endif else begin
            ; blank filter selected
            device, set_graphics = 6
            thick = 3
            col = off_col
            if ( i eq curbox ) then col = sel_col
            xyouts, box_x(i), box_y(i)-15, strcompress(i, $
              /remove_all), charthick=thick, color=col, charsize=4, $
              /device, alignment=0.5
            device, set_graphics = 3
        endelse
    endif
endif
endfor


END
;**********************************************************

;**********************************************************
PRO erase_curbox
common shared_1
common shared_2
common shared_3
common shared_4
common shared_5
common shared_6
common shared_ao

widget_control, draw, get_value=index
wset, index
if ( boxes(curbox).spmode eq 'high res' ) then begin
    device, set_graphics = 6
    thick = 3
    col = sel_col
    plots, [boxes(curbox).x(0),boxes(curbox).x(1)], $
      [boxes(curbox).y(0),boxes(curbox).y(1)], $
      color = col, thick = thick
    plots, [boxes(curbox).x(2),boxes(curbox).x(1)], $
      [boxes(curbox).y(2),boxes(curbox).y(1)], $
      color = col, thick = thick
    plots, [boxes(curbox).x(2),boxes(curbox).x(3)], $
      [boxes(curbox).y(2),boxes(curbox).y(3)], $
      color = col, thick = thick
    plots, [boxes(curbox).x(0),boxes(curbox).x(3)], $
      [boxes(curbox).y(0),boxes(curbox).y(3)], $
      color = col, thick = thick
    boxlab = strtrim(string(curbox),1)
    xlab = 0.1*boxes(curbox).x(1) + 0.9*boxes(curbox).x(0)
    ylab = 0.9*boxes(curbox).y(0) + 0.1*boxes(curbox).y(3)
    xyouts, xlab, ylab, boxlab, $
      charsize = 4, color = col, charthick = 3
    device, set_graphics = 3
endif
if (boxes(curbox).spmode eq 'low res') then BEGIN
    device, set_graphics = 6
    thick = 3
    col = sel_col
    plots, [boxes(curbox).x(0),boxes(curbox).x(1)], $
      [boxes(curbox).y(0),boxes(curbox).y(1)], $
      color = col, thick = thick
    plots, [boxes(curbox).x(2),boxes(curbox).x(1)], $
      [boxes(curbox).y(2),boxes(curbox).y(1)], $
      color = col, thick = thick
    plots, [boxes(curbox).x(2),boxes(curbox).x(3)], $
      [boxes(curbox).y(2),boxes(curbox).y(3)], $
      color = col, thick = thick
    plots, [boxes(curbox).x(0),boxes(curbox).x(3)], $
      [boxes(curbox).y(0),boxes(curbox).y(3)], $
      color = col, thick = thick
    boxlab = strtrim(string(curbox),1)
    xlab = 0.1*boxes(curbox).x(1) + 0.9*boxes(curbox).x(0)
    ylab = 0.9*boxes(curbox).y(0) + 0.1*boxes(curbox).y(3)
    xyouts, xlab, ylab, boxlab, $
      charsize = 4, color = col, charthick = 3
    device, set_graphics = 3
ENDif
if boxes(curbox).spmode eq 'imaging' then BEGIN
;    plot, [1,46], [1,46], xtitle='arcseconds', $
;      ytitle='arcseconds', xstyle=1, ystyle=1, /nodata
;    xyouts, 23, 23, 'IMAGING', /data, charsize=8, align=0.5, $
;      charthick=10
;    boxes(curbox).spmode = 'imaging'
    device, set_graphics = 6
    thick = 3
    col = sel_col
    xyouts, box_x(curbox), box_y(curbox)-15, strcompress(curbox, $
      /remove_all), charthick=thick, color=col, charsize=4, $
      /device, alignment=0.5
    device, set_graphics = 3
ENDif

END
;**********************************************************

;**********************************************************
PRO format_quit
; Destroys window when quit is selected
common shared_1
common shared_2
common shared_3
common shared_4
common shared_5
common shared_6
common shared_ao

widget_control, formatbase, /destroy
END

;**********************************************************
;             Begin User Lines Procedures
;** *** *** *** *** *** *** *** *** *** *** *** *** *** ***

pro browse_button_event, event

common shared_1
common shared_2
common shared_3
common shared_4
common shared_5
common shared_6	
common shared_ao

        ; if list has changed, changed=1

        ; if list has changed
        if changed then begin
                message=['List contents have changed.', $
			'Do you wish to save list first?']
		answer=dialog_message(message, dialog_parent=event.id, $
			/question)
		if answer eq 'Yes' then begin
			save_as_button_event, event
		endif
	endif

        ; open pick file dialog
	filelist=dialog_pickfile(file=current_dir+'/'+'user.lst', $ 
		filter='*.lst', /read, group=event.id)
	if filelist ne '' then begin
		open_list
	endif
end

pro add_line_button_event, event

common shared_1
common shared_2
common shared_3
common shared_4
common shared_5
common shared_6
common shared_ao


	widget_control, add_line_box, get_value=new_line
	widget_control, add_comment_box, get_value=new_comment

	if new_line[0] ne '' then begin   ; make sure line isn't empty
		if list[0] ne '' then begin  ; if list isn't already empty
			list_size=size(list)
			if list_size(0) eq 0 then begin ;if just one line
				list=[list, new_line]
                                comments=[comments, new_comment]
				order=sort(list)
                                list=list(order)
                                comments=comments(order)
				displayed_list=string(list)+'     '+comments
			endif else begin  ; if more than one line
				temp_array=list
				temp_comments=comments
				list=fltarr(list_size(1)+1)
				comments=strarr(list_size(1)+1)
				list(list_size(1))=new_line
				comments(list_size(1))=new_comment
				list[0:list_size(1)-1]=temp_array
				comments[0:list_size(1)-1]=temp_comments
				order=sort(list)
				list=list(order)
				comments=comments(order)
				displayed_list=string(list)+'     '+comments
			endelse
		endif else begin     ; if empty list
			list=new_line
			comments=new_comment
			displayed_list=string(list)+'     '+comments
		endelse
	endif

        changed=1
	widget_control, line_list, set_value=displayed_list

end

pro rgb_event, event

common shared_1
common shared_2
common shared_3
common shared_4
common shared_5
common shared_6
common shared_ao
	
	base_id = widget_info(event.id, /parent)   ;get slider id's

        ; get new value of the slider
	if event.id eq base_id+1 then r=event.value
	if event.id eq base_id+2 then g=event.value
	if event.id eq base_id+3 then b=event.value

	user_col=r+256L*(g+256L*b)

        wset, color_id
	polyfill, [0,1,1,0], [0,0,1,1], /line_fill, /normal, $  ;fill box
		color=user_col
end

pro ok_button_event, event

common shared_1
common shared_2
common shared_3
common shared_4
common shared_5
common shared_6
common shared_ao

        ; update box on user lines dialog
        user_col=r+256L*(g+256L*b)

        wset, draw_id
	polyfill, [0,1,1,0], [0,0,1,1], /line_fill, /normal, color=user_col
	
        widget_control, color_base, /destroy
        color_base=0
end
	

pro remove_buttons_event, event

common shared_1
common shared_2
common shared_3
common shared_4
common shared_5
common shared_6
common shared_ao

	case event.value of
		0: begin  ; remove line button
                        changed=1
                        ; get selection and size of selection
                        ; on unix, selection size is limited to one
			info=widget_info(line_list, /list_select) 
			;print, info, info[0]
			info_size=size(info)
			;print, 'Info size=', info_size
			temp_array=list
			temp_comments=comments
			list_size=size(list)
                        if list_size(0) eq 0 then i_max=0 $
                          else i_max=list_size(1)-1
			;print, list_size
			if info[0] ge 0 then begin
				if info_size[0] eq 0 then begin
					j=0
					for i=0, i_max do begin
						if i ne info[0] then begin
							temp_array[j]=list[i]
						temp_comments[j]=comments[i]
							j=j+1
						endif
					endfor
					if j gt 0 then begin
						list=fltarr(j)
						comments=strarr(j)
						list[*]=temp_array[0:j-1]
					comments[*]=temp_comments[0:j-1]
					endif else begin
						list=''
						comments=''
					endelse
				endif else begin
					n=0
					m=0
					for l=0, list_size(1)-1 do begin
						if l ne info[m] then begin
							temp_array[n]=list[l]
						temp_comments[n]=comments[l]
							n=n+1
						endif else begin
						       m=(m+1<(info_size(1)-1))
						endelse
					endfor
					if n gt 0 then begin
						list=fltarr(n)
						comments=strarr(n)
						list[*]=temp_array[0:n-1]
					comments[*]=temp_comments[0:n-1]
					endif else begin
						list=''
						comments=''
					endelse
				endelse
			endif
			displayed_list=string(list)+'     '+comments
			widget_control, line_list, set_value=displayed_list
		end

		1: begin  ; clear list button
                        changed=1
			warning=["This will erase all","lines from your list!"]
			answer=dialog_message(warning, /default_cancel, $
				dialog_parent=event.id, /cancel)
			if answer eq 'OK' then begin
				list=''
				comments=''
				widget_control, line_list, set_value=list
			endif
		end

		2: begin  ; select color button
			
                        if color_base ne 0 then begin
                            widget_control, color_base, /show
                        endif else begin
                            ; create new widget base and controls
                            color_base=widget_base(/col, group_leader= $
                                user_base, title='Select Color')
                            color_sel=widget_draw(color_base, xs=275, ys=30)
                            rgb=cw_rgbslider(color_base)
                            ok_button=widget_button(color_base, value='OK', $
                                     font=font)
                            widget_control, color_base, /realize
                            widget_control, color_sel, get_value=color_id
                            ; set slider initial values
                            widget_control, rgb+10, set_value=r
                            widget_control, rgb+11, set_value=g
                            widget_control, rgb+12, set_value=b
                            
                            ;fill box on dialog with initial color
                            wset, color_id
                            polyfill, [0,1,1,0], [0,0,1,1], /line_fill, $
                                  /normal, color=user_col
                            
                            xmanager, 'ok_button', ok_button, /just_reg
                            xmanager, 'rgb', rgb, /just_reg
                            xmanager

                        endelse
		end
	endcase

end


pro save_as_button_event, event

common shared_1
common shared_2
common shared_3
common shared_4
common shared_5
common shared_6
common shared_ao

	widget_control, save_list_box, get_value=save_list_name
        
        ; check to see is edit box begins with a slash i.e. a directory
	if strmid(save_list_name[0], 0, 1) eq '/' then begin
		default=save_list_name[0] 
	endif else begin  ; if not, add current path
		default=current_dir+'/'+save_list_name[0]
	endelse

        ; bring up save as dialog
	save_file=dialog_pickfile(/write,  file=default, filter='*.lst', $
		group=event.id)
	print, save_file
	if save_file ne '' then begin
		result=findfile(save_file)  ; see if file already exists
		if result[0] eq '' then begin
			save_list
		endif else begin  ; if yes, prompt for overwrite
			message=['File ' + save_file+' exists.', $
				'Do you wish to Overwrite?']
			answer=dialog_message(message,title='Overwrite File?',$
				/default_no, dialog_parent=event.id, /question)
			if answer eq 'Yes' then save_list
		endelse
	endif
end

pro control_buttons_event, event

common shared_1
common shared_2
common shared_3
common shared_4
common shared_5
common shared_6
common shared_ao


        ; get list size
        user_list_size=size(list)
        if user_list_size(0) eq 0 then begin
                user_list_length=1
        endif else begin 
                user_list_length=user_list_size(1)
        endelse

        ; get value of z and update redshift value
        widget_control, z_box, get_value=redshift
        ; create new list array to be plotted on echellogram
        shifted_list=(1+float(redshift(0)))*list
        ; print, shifted_list
        case event.value of
                0: begin   ; ok button
                    if changed then begin
                          message=['List contents have changed.', $
				'Do you wish to save list first?']
			answer=dialog_message(message, dialog_parent=event.id,$
				/question)
			if answer eq 'Yes' then begin
				save_as_button_event, event
			endif
                    endif
                    draw_echellogram
                    draw_boxes
                    widget_control, event.top, /destroy
                    user_base=0
                end

		1:begin   ; apply button
                    draw_echellogram
                    draw_boxes
                end

                2: begin  ; cancel button
                    list=old_list
                    comments=old_comments
                    widget_control, user_lines_id, set_value='[  ]   User Lines  <Ctrl+u>'
                    widget_control, i_menu+6, set_value='user'
                    user_lines_setup
                    widget_control, event.top, /destroy
                    user_base=0
		end

                else: print, event.value
	endcase
end

pro open_list

common shared_1
common shared_2
common shared_3
common shared_4
common shared_5
common shared_6
common shared_ao

        ; create temp array since we don't know the size of data list
	flt_array=fltarr(1000)
	temp_comments=strarr(1000)

        ; call read_lines procedure
	read_lines, filelist, flt_array, temp_comments, i, err 

        ; if file does not exist
	if err eq -215 then begin 
		message=['File does not exist', $
			 'Do you wish to create a new one?']
		answer=dialog_message(message, dialog_parent=browse_button, $
			/question)
		if answer eq 'Yes' then begin
			list=''
                        comments=''
                        changed=0
			widget_control, filelist_box, set_value=filelist
			widget_control, save_list_box, set_value=filelist
		endif
	endif else begin
            changed=0
            
            ; extract correctly sized list from temp arrays
            if i ge 1 then begin
                list=fltarr(i)
		comments=strarr(i)
		for j=0, i-1 do begin
			list(j)=flt_array(j)
			comments(j)=temp_comments(j)
		endfor
            endif else begin    ;if temp arrays are empty, so are the new lists
                list=''
                comments=''
            endelse
            
            widget_control, filelist_box, set_value=filelist
            widget_control, save_list_box, set_value=filelist
            
            ; update the old lists for checking later
            old_list=list
            old_comments=comments

            displayed_list=string(list)+'     '+comments
            widget_control, line_list, set_value=displayed_list
        endelse
end

pro save_list

common shared_1
common shared_2
common shared_3
common shared_4
common shared_5
common shared_6
common shared_ao


	list_size=size(list)
	openw, 1, save_file, error=err
	print, err

        ;if error saving file (directory does not exist)
	if err eq -215 then begin
		message=['ERROR: Cannot save file.', 'Check filename.']
		answer=dialog_message(message, dialog_parent=save_list_box)
;		if answer eq 'Yes' then begin
;			sep=str_sep(save_list_name[0], '/')
;			sep_size=size(sep)
;			print, sep_size
;				if sep_size(1) eq 2 then begin
;					directory=sep(0)
;					command='mkdir ' + directory
;					spawn, command
;					new_message=['Directory '+directory+ $
;					     '/ created.', $;
;					     'Press save again to save file.'] 
;					new=dialog_message(new_message, $
;						dialog_parent=save_list_box, $
;						/information)
;				endif else begin
;				new_message=['Cannot create directory.']
;				new=dialog_message(new_message, /information, $
;					dialog_parent=save_list_box)
;				endelse
;		endif
	endif else begin  ; if no error
		if list[0] eq '' then begin  ; if empty list
			printf, 1, ''
		endif else begin 
			if list_size[0] eq 0 then begin   ; if one item in list
				printf, 1, list+' '+comments
			endif else begin  ; if multiple items in list
				for i=0, list_size(1)-1 do begin
					line=string(list(i))+' '+comments(i)
					printf, 1, line
				endfor
			endelse
		endelse
		close, 1

                ; update old lists
		old_list=list
		old_comments=comments

                changed=0

                widget_control, filelist_box, set_value=save_file
                widget_control, save_list_box, set_value=save_file
            endelse

end

pro user_lines

common shared_1
common shared_2
common shared_3
common shared_4
common shared_5
common shared_6

print, 'in user_lines.pro..., list=',list, ', color=', user_col

changed=0
color_base=0
remove=['Remove Selection', 'Clear List', 'Select Color']
control=['OK', 'Apply', 'Cancel']


if user_base ne 0 then begin
    widget_control, user_base, /show
endif else begin

; create widgets
user_base=widget_base(/col, Title='User Line List')
file_base=widget_base(user_base, /row)
filelist_box=cw_field(file_base, value=filelist, title='List Filename:', $
	xs=30, /noedit, font=font)
browse_button=widget_button(file_base, value='Browse', font=font)
line_list=widget_list(user_base, value=displayed_list, ys=4)
list_text=widget_label(user_base, value='List line values in microns.')
add_line_base=widget_base(user_base, /row)
add_line_box=cw_field(add_line_base, value=new_line, /floating, $
	title='Add Line to List:', xs=30, font=font)
add_line_button=widget_button(add_line_base, value='Add', font=font)
add_comment_box=cw_field(user_base, value='', title='New Line Comment:', xs=35, $
                        font=font)
z_box=cw_field(user_base, value=redshift, title='Z = ', font=font)
remove_base=widget_base(user_base, /row)
remove_buttons=cw_bgroup(remove_base, remove, /row, font=font)
draw_base=widget_base(remove_base, /col)
color_draw=widget_draw(draw_base, xs=60, ys=20)
save_base=widget_base(user_base, /row)
save_list_box=cw_field(save_base, value='No File Opened...', xs=30, $
	title='Save List As:', font=font) 
save_as_button=widget_button(save_base, value='Save', font=font)
control_buttons=cw_bgroup(user_base, control, /row, font=font)

widget_control, user_base, /realize
widget_control, color_draw, get_value=draw_id

;fill box with initial color
wset, draw_id
polyfill, [0,1,1,0], [0,0,1,1], /line_fill, /normal, color=user_col

xmanager, 'browse_button', browse_button, /just_reg, /no_block
xmanager, 'add_line_button', add_line_button, /just_reg, /no_block
xmanager, 'remove_buttons', remove_buttons, /just_reg, /no_block
xmanager, 'save_as_button', save_as_button, /just_reg, /no_block
xmanager, 'control_buttons', control_buttons, /just_reg, /no_block
xmanager

endelse

end

;** *** *** *** *** *** *** *** *** *** *** *** *** *** ***
;               End User Lines Procedures
;**********************************************************


PRO nirspec_efs, server=server, ao=ao

common shared_1
common shared_2
common shared_3
common shared_4
common shared_5
common shared_6
common shared_ao
common shared_go
common plot_shared

if keyword_set (ao) then begin
 aomode=1
endif else begin
 aomode=0
endelse

if not keyword_set(server) then server='none'
case server of 
    'nirspec': serv='nirspec'
    'nirspecsim': serv='nirspecsim'
    'sim': serv='nirspecsim'
    else: serv='none'
endcase


!except=2
;----------------------------
; setup defaults
;----------------------------
font = '-adobe-times-bold-r-normal--14-100-100-100-p-76-iso8859-1'
;    font = '-adobe-helvetica-bold-o-normal--20-140-100-100-p-103-iso8859-1'
;    font = '-adobe-times-bold-r-normal--20-140-100-100-p-100-iso8859-1'
num_scr = 0			; script number = 0

date = string(systime())        ; Get the system date.
dd = strmid(date, 8, 2)
mm = strmid(date, 4, 3)

scr_name_base = 's' +  string(format='(I2.2)',dd) + mm
scriptname = scr_name_base + string(format='(I4.4)',num_scr) + '.csh'
;scriptdir = '/kroot/kss/nirspec/ui/xnirspec/scripts'
;spawn, 'show -s nirspec -terse scriptdir', scriptdir
scriptdir2 = '.'

nodmode       = 1               ; nod 1 stare
calmode       = 0               ; setup only
plotted=0                       ; flag for if filter plot window is up.
blank=0                         ; flag if blank is selected filter
lock_flag=0
dlambda = fltarr(1000)          ; wavelengths for sampling
echphi = fltarr(1000)
graphi = fltarr(1000)
lgraphi = fltarr(1000)
hgraphi = fltarr(1000)
leftlam = fltarr(100)
rightlam = fltarr(100)
leftech = fltarr(100)
rightech = fltarr(100)
lleftx = fltarr(100)
lrightx = fltarr(100)
llefty = fltarr(100)
lrighty = fltarr(100)
uleftx = fltarr(100)
urightx = fltarr(100)
ulefty = fltarr(100)
urighty = fltarr(100)
leftgra = fltarr(100)
rightgra = fltarr(100)
numorders = fix(10)
;-----------------------------
; Setup wavelength overlays
neon_stat = 0                   ; Turn off neon lines
krypton_stat = 0                ; Turn off krypton lines
argon_stat = 0                  ; Turn off argon lines
xenon_stat = 0                  ; Turn off xenon lines
private_stat = 0                ; Turn off private line list
oh_stat = 0                     ; Turn off OH lines.
user_stat = 0                   ; Turn off user lines.

neon_num = 70
neon_lines = fltarr(neon_num)
argon_num = 124
argon_lines = fltarr(argon_num)
krypton_num = 91
krypton_lines = fltarr(krypton_num)
xenon_num = 57
xenon_lines = fltarr(xenon_num)
oh_num = 65
oh_lines = fltarr(oh_num)
user_list_length=0

;-----------------------------
; User lines initial values
user_base=0
filelist='Click browse to open file...'
;list=[0.0008, 0.0009, .00102, .00105, .00114, .00115, .00118]
;comments=strarr(7)
list=''
comments=''
displayed_list=string(list)+comments
old_list=list
old_comments=comments
new_line=''
cd, '.', current=current_dir
r=255 & g=0 & b=255 ;inital color is purple
result=0
redshift=0
working_dir = '~/setups'

;----------------------------
; Include the ns_params.pro file to set the gratings and detector angles
;----------------------------

if (aomode) then begin
 @ns_params_ao.pro
endif else begin
 @ns_params.pro
endelse

; @/kroot/kss/nirspec/ui/idllib/ns_params.pro
; Set up echelle parameters
echtheta = 63.00*!PI/180.       ; Echelle angle
;echsigma = 1./(23.21*1.00396)     ; echelle groove density
;echsigma = 1./(23.0)           ; echelle groove density: match spectra 3/6/98
;echsigma = 1./(22.74)          ; echelle groove density: match spectra 3/6/98
;echgamma = (5./180.0)*!PI       ; the QLM angle
;---------------------------------
; Some cross disperser parameters
;grasigma = 1./(75.465*1.00396)      ; cross disperser groove density
gratheta = 10.0*!pi/180.0        ; cross disperser blaze
;grahaang = 24.96*!pi/180.        ; the half opening angle
;grahaang = 25.00*!pi/180.        ; the half opening angle
;--------------------------------
; The TMA 
;few = 465                       ; effective focal length in slit width
;feh = 406                       ; effective focal length in slit height
;dist_few = fltarr(4)
;dist_feh = fltarr(4)
;dist_few = [491, 447, 442, 492] ; effective fl in slit wid. for corners
;dist_feh = [425, 405, 399, 425] ; effective fl in slit wid. for corners
;--------------------------------
; The Detector
;pix = 0.027			; pixel size in mm
;npix=1024			; size of array
;yarc = 0.2			; arcseconds per y pixel in high or lowres
;xarc = 0.16			; arcseconds per x pixel in high res
;chipang = 4.657                 ; Rotation angle of chip.

;----------------------------
; Set-up Parameters
moving = 0                      ; moving =1 if user is dragging box.
maxbox = 1000                   ; Maximum number of grating positions allowed
curbox = 1                      ; box currently being defined
totbox = 1                      ; Number of boxes currently defined
A = {position, $                ; For each position (box) setup structure of param.
     echphi: 0.0, $             ; Echelle Angle
     graphi: 0.0, $             ; Cross disperser angle
     spmode: 'high res', $      ; observing mode
     reps: 1, $                 ; Number of reps per nod pattern
     slit: 0, $                 ; Slit position
     filter: filterarr(0).position, $	; Filter position
     blocker: 1, $              ; Blocking filter
     wheel: filterarr(0).wheel, $	; Which filter wheel name
     object: ' ', $             ; Object name
     otint: 10.0, $             ; Integration time on object
     star: ' ', $               ; Stars name
     stint: 5.0, $              ; Integration time on star
     rot: 45.0, $               ; Rotator position for all exposures
     coadds: 1, $               ; Number of Coadds
     scoadds: 1, $              ; Number of Coadds
     x: fltarr(4), $            ; Array containing x coords of box corners
     y: fltarr(4) $             ; Array containing y coords of box corners
    }
boxes = REPLICATE({position},maxbox) ; Memory for all positions	
; Setup initial positions of 1st position
boxes(curbox) = {position, echtheta, gratheta, 'high res', 1, 0, $
                 filterarr(0).position, 1, filterarr(0).wheel, $
                 'default', 10.0, 'default', 5.0, 92.0, 1, 1, fltarr(4), $
                 fltarr(4)}
boxes(0)=boxes(curbox)
slitlen = slitlist(boxes(curbox).slit).length

filtermenu_index = intarr(maxbox)   ; index of the filter menu

; get bad pixel mask...
mask_flag=0
; read in fits file
;bad_array=readfits('mask2.fits', hd)
;bad_array=readfits('nirspec_bad.fits', hd)
;bad_array=readfits('bad.fits', hd)

temp_bad_x=intarr(12000)
temp_bad_y=intarr(12000)

read_data, 'bad.lst', temp_bad_x, temp_bad_y, numbad

; check to make sure file is valid

; get number of bad pixels.  since each bad pixel has a value of one, 
; and valid pixels have a value of zero, the number of bad pixels 
; is equal to the sum of the values of the array.
numbad=numbad>1
; define structure to hold bad pixels
bad={pixel, $
     x:0.0, $
     y:0.0, $
     echdel:0.0, $
     gradel:0.0, $
     plot_x:0.0, $
     plot_y:0.0}
bad=replicate({pixel}, numbad)
; extract bad pixels and fill structure
;c=0
;for i=0, 1023 do begin
;    for j=0, 1023 do begin
;        if bad_array(i, j) then begin
;            bad(c).x=i
;            bad(c).y=j
;            c=c+1
;        endif
;    endfor
;endfor

bad.x=temp_bad_x[0:numbad-1]
bad.y=temp_bad_y[0:numbad-1]

; for each pixel, find rotated position and echelle and cross
; dispersor angles
chiprot=!pi*chipang/180.0
st=sin(chiprot)
ct=cos(chiprot)
for i=0, numbad-1 do begin
    bad(i).x=bad(i).x-512
    bad(i).y=bad(i).y-512
    ; the following arguments are such to account for rotation of echellogram
    bad_few=neflw(bad(i).y, -bad(i).x)
    bad_feh=neflh(bad(i).y, -bad(i).x)
    tx=bad(i).x*ct-bad(i).y*st
    ty=bad(i).x*st+bad(i).y*ct
    bad(i).x=tx
    bad(i).y=ty
    bad(i).echdel=atan(bad(i).x*pix/bad_few)
    bad(i).gradel=atan(bad(i).y*pix/bad_feh)
endfor

;loadct, 39                      ; load Rainbow+white table
;loadct, 3                      ; load red heat table
;!P.NOCLIP=1
print, "number of colors =", !d.n_colors

;----------------------------
; SETUP PULLDOWN MENUS
;----------------------------
junk = { CW_PDMENU_S, flags:0, name:'' }
format_menu_desc = [ { CW_PDMENU_S, 1, 'FILE' }, $
                     { CW_PDMENU_S, 0, 'Open Configuration...     <Ctrl+o>' },$
                     { CW_PDMENU_S, 0, 'Save Configuration As... <Ctrl+s>' }, $
		     { CW_PDMENU_S, 0, 'Configure Script File       <Ctrl+f>'}, $
                     { CW_PDMENU_S, 2, $
                       'Quit                                       <Ctrl+q>' }, $ 
                     { CW_PDMENU_S, 1, 'TELESCOPE' }, $
;                     { CW_PDMENU_S, 0, 'XLoadCT...' }, $
;                     { CW_PDMENU_S, 2, 'Set Label Offsets...' }, $
                     { CW_PDMENU_S, 2, 'Set Nod Parameters  <Ctrl+p>' }, $
                     { CW_PDMENU_S, 1, 'OVERLAYS' }, $
                     { CW_PDMENU_S, 0, 'Clear Overlays  <Ctrl+c>' }, $
                     { CW_PDMENU_S, 0, '[  ]   User Lines  <Ctrl+u>' }, $
;                     { CW_PDMENU_S, 0, 'Wavelength...' }, $
                     { CW_PDMENU_S, 0, '[  ]   OH Lines    <Ctrl+l>' }, $
;                     { CW_PDMENU_S, 0, 'Edge lambda'}, $
                     { CW_PDMENU_S, 3, 'Arc Lines'}, $
                     { CW_PDMENU_S, 0, '[  ]   Neon         <Ctrl+n>'}, $
                     { CW_PDMENU_S, 0, '[  ]   Argon       <Ctrl+a>'}, $
                     { CW_PDMENU_S, 0, '[  ]   Krypton  <Ctrl+k>'}, $
                     { CW_PDMENU_S, 2, '[  ]   Xenon       <Ctrl+x>'}, $
                     { CW_PDMENU_S, 3, 'HELP' }, $
                     { CW_PDMENU_S, 0, 'Keyboard Hotkey List'}, $
                     { CW_PDMENU_S, 2, 'User Manual' } ]

;----------------------------
; Layout widget bases
;----------------------------
if (aomode) then begin
formatbase = widget_base(/col, mbar=base_mbar, /TLB_SIZE_EVENTS, $
                         TITLE='NIRSPEC Format Simulator (AO Mode)', $
                         resource_name='efsbase')
endif else begin
formatbase = widget_base(/col, mbar=base_mbar, /TLB_SIZE_EVENTS, $
                         TITLE='NIRSPEC Format Simulator', $
                         resource_name='efsbase')
endelse
formatmenu = CW_PDMENU(base_mbar, format_menu_desc, /RETURN_NAME, /MBAR, $
                       FONT=font)
upperbase = widget_base(formatbase, /row)
leftbase = widget_base(upperbase, /col)
drawbase = widget_base(upperbase, /col)
rightbase = widget_base(upperbase, /col)
boxbase = widget_base(drawbase, /row)

; copy wid of format base to be passed to plot filter
efs_base=formatbase

; set initial screen sizes
;dx = 512
;dy = 512
scrx = 350
scry = 350
plot_x1 = 30
plot_x2 = scrx-30
plot_y1 = 30
plot_y2 = scry-30

; jlw box moving inits
sptype='high res'               ; set by specmenu
current_fil=0                   ; set by filtermenu
current_slit=0                  ; set by slitmenu
valid_boxes=intarr(maxbox)      ; array for boxes containing point clicked
csr_x=175                       ; cursor x location
csr_y=175                       ; cursor y location
in_draw=0                       ; flag if cursor is in draw window
xoff=0                          ; x offset of click from center of box
yoff=0                          ; y offset of click from center of box
box_x=intarr(maxbox)            ; x locs of center of boxes (device) 
box_y=intarr(maxbox)            ; y locs of center of boxes (device)
box_x(curbox)=175
box_y(curbox)=175

; Setup the colors for the display.
;RED =   [0, 255,  50,  50, 200, 255]
;GREEN = [0, 160, 255,  50,  50, 255]
;BLUE =  [0,   0,   0, 255,   0, 255]
;tvlct, RED, GREEN, BLUE
;sel_col = 1
;off_col = 2
;spec_col = 3
;border_col = 4
;text_col = 5
;if ( !d.n_colors ge 16000000 ) then begin		; 24 bit color
	; colors are R + 256L*(G+ 256L*B)
	sel_col = 255L + 256L*(160L+256L*0L)	; Orange
	off_col = 50L + 256L*(255L+256L*0L)	; Green
	spec_col = 50L + 256L*(50L+256L*255L)	; Blue
	border_col = 200L + 256L*(50L + 256L*0L) ; Red
	text_col = 255L + 256L*(255L+256L*255L)	; White
        neon_col = 255L + 256L*(200L+256L*50L)  ; Orange
        argon_col = 220L + 256L*(220L+256L*220L); Gray
        krypton_col = 255L + 256L*(100L+256L*0L) ; Blue
        xenon_col = 100L + 256L*(255L+256L*100L) ; Green
        oh_col = 200L + 256L*(200L+256L*0L)     ; Yellow
        user_col = 255L + 256L*(0L+256L*255L)   ; Purple
;	print, '24 bit mode'
;endif
;if ( !d.n_colors lt 16000000 ) then begin
;	sel_col = 80.0/100.0 * !d.n_colors
;	off_col = 36.0/100.0 * !d.n_colors
;	spec_col = 25.0/100.0 * !d.n_colors
;	border_col = 95.0/100.0 * !d.n_colors
;	text_col = !d.n_colors - 1
;	print, "less that 24bits"
;endif

draw = widget_draw(drawbase, xsize=scrx, ysize=scry, /button_events, $
                  /motion_events, /tracking_events)

d2base = widget_base(drawbase,/ROW)
lamlab = widget_label(d2base, value='Lambda:   ', font=font)
lambox = widget_label(d2base, value='0', font=font, xsize=80, /align_right, $
                     frame=2)
ordlab = widget_label(d2base, value='     Ech Order: ', font=font)
ordbox = widget_label(d2base, value='0', font=font, xsize=80, /align_right, $
                     frame=2)

d3base = widget_base(drawbase,/ROW)
xpixlab = widget_label(d3base, value='X-PIXEL:', font=font)
xpixbox = widget_label(d3base, value='0', font=font, xsize=80, /align_right, $
                     frame=2)
ypixlab = widget_label(d3base, value='     Y-PIXEL:  ', font=font)
ypixbox = widget_label(d3base, value='0', font=font, xsize=80, /align_right, $
                     frame=2)

;d1base = widget_base(drawbase,/ROW)
;echlab = widget_label(d1base, value='Echelle:    ', font=font)
;echbox = widget_label(d1base, value='0', font=font, xsize=80, /align_right, $
;                     frame=2)
;echbox = cw_field(d1base, TITLE='Echelle:', FONT=font,/row, $
;                 xs=10, value = 180.0*boxes(curbox).echphi/!pi)

;gralab = widget_label(d1base, value='     Cross Disp:', font=font)
;grabox = widget_label(d1base, value='0', font=font, xsize=80, /align_right, $
;                     frame=2)
;grabox = cw_field(d1base, TITLE='Cross Disp:', FONT=font,/row, $
;                 xs=10, value = 180.0*(boxes(curbox).graphi+grahaang)/!pi)


fullname = ['Setup Only', 'Lamps Only', 'Object Only','Bright Object Only'] 
; 'Object+Lamps', 'Object+Star', 'Full Sequence']
full = WIDGET_DROPLIST(leftbase, value = fullname, title='Obs. Mode', $
                       font = font, /align_left)

nodname = ['Stare', 'Nod 2', 'Nod 3', 'Nod 4', 'ABBA', 'Nod Off', 'User']
nods = WIDGET_DROPLIST(leftbase, value = nodname, title='Nod Pattern', $
                       font = font, /align_left)

repsname = ['One', 'Two', 'Three', 'Four', 'Five', 'Six']
reps = WIDGET_DROPLIST(leftbase, value = repsname, title='Number of Reps', $
                       font = font, /align_left)

specname = ['High Res', 'Low Res', 'Imaging']
specmode = WIDGET_DROPLIST(leftbase, value = specname, title='Spec. Mode', $
                           font = font, /align_left)

slitmenu = WIDGET_DROPLIST(leftbase, value = slitlist.name, title='Slit', $
                           font = font, /align_left)

filtermenu = WIDGET_DROPLIST(leftbase, value = filterarr.name, title='Filter',$
                             font = font, /align_left)

;blocking_text=widget_label(leftbase, value='Blocking:', font=font, $
;                          /align_left)
if not (aomode) then begin
 thinname=['Open', 'Thin']
 thinmenu=cw_toggle(leftbase, thinname, title='Blocking:', xs=150, font=font, $
                   value=1, event_pro='thinmenu_event')
endif 
;thinmenu= cw_bgroup(leftbase, thinname, /exclusive, set_value=1, font=font, $
;                    row=1, /frame)
; see notes about cw_toggle below

plot_filter_button=widget_button(leftbase, value='Plot Filter Profile', $
                                  font=font, /align_center, xs=200)

; for cw_toggle, make sure event_pro keyword is set to indicated
; handler.  Also do not add this button to xmanager.  Also, note 
; if-loop in formatbase_event. since this event is not in xmanager
; the formatbase event is called whenever an event occurs with this 
; button.  add the loop so that resize events only occur when the
; event acts on the top level base.  otherwise, do nothing and
; call event_pro indicated below.
mask_toggle=cw_toggle(leftbase, ['Show', 'Hide'], title='Bad Pixel Mask:', $
                      font=font, event_pro='mask_toggle_event', value=1)

echbox = cw_field(leftbase, TITLE='Echelle:', FONT=font,/row, $
                 xs=10, value = round2(180.0*boxes(curbox).echphi/!pi))
grabox = cw_field(leftbase, TITLE='Cross Disp:', FONT=font,/row, $
                 xs=10, value = round2(180.0*(boxes(curbox).graphi+grahaang)/!pi))

lock_toggle=cw_toggle(leftbase, ['Unlock', 'Lock'], title='Lock Mechanisms:', $
                      font=font, event_pro='lock_toggle_event', value=0)



;curbox_base=widget_base(leftbase, /row)
;curbox_button=widget_button(curbox_base, value='Current Box:', font=font)
;curbox_field=widget_text(curbox_base, value='1', xs=4)



;set_label = widget_label(leftbase, font=font, $
;                         value=' ')

box_base=widget_base(formatbase, /row)
set_label = widget_label(box_base, font=font, $
                         value='SET-UP BOXES:')
box_button_base=widget_base(box_base, /row, frame=2)
; button size changed by JLW 2/24/99, old val=110
new_box = widget_button(box_button_base, val='NEW', font=font, xs=85, $
                        /align_center)
pre_box = widget_button(box_button_base, val='PREV', font=font, xs=85, $
                        /align_center)
next_box = widget_button(box_button_base, val='NEXT', font=font, xs=85, $
                         /align_center)
del_box = widget_button(box_button_base, val='DELETE', font=font, xs=85, $
                        /align_center)
curbox_button=widget_button(box_base, value='Current:', font=font)
curbox_field=widget_text(box_base, value=strcompress(curbox), xs=4)

print_preview = 0
;----------------------------------
; build the non-script user interface
;----------------------------------
ui_outer_base = widget_base(formatbase, /COL, frame=4)
;uilabel = widget_label(ui_outer_base, font=font, $
;	value='EXPOSURE INFORMATION')
ui_base = widget_base(ui_outer_base, /ROW)

name_base = widget_base(ui_base, /COL)
obj_name = cw_field(name_base, TITLE='Object', FONT=font,/ROW, $
                    value = boxes(curbox).object, xs=20)
star_name = cw_field(name_base, TITLE='Star', FONT=font,/ROW,$
                     value = boxes(curbox).star, xs=20)
itime_base = widget_base(ui_base, /col)
itime = cw_field(itime_base, TITLE='Itime', FONT=font,/row, $
                 xs=10, value = boxes(curbox).otint)
stime = cw_field(itime_base, TITLE='Itime', FONT=font,/row, $
                 xs=10, value = boxes(curbox).stint)

coadds_base = widget_base(ui_base, /col)
coadds = cw_field(coadds_base, TITLE='Coadds', FONT=font,/row,xs=3, $
                  value = boxes(curbox).coadds)
scoadds = cw_field(coadds_base, TITLE='Coadds', FONT=font,/row,xs=3, $
                   value = boxes(curbox).scoadds)
;boxnum = cw_field(ui_base, TITLE='Box', val=curbox, FONT=font, /row, xs=3)

go_button = widget_button(ui_base, VALUE='   GO   ', $
                          resource_name='test_go', FONT=font, ys=75, xs=60)

abort_button = widget_button(ui_base, VALUE='ABORT', $
                             resource_name='test_abort', FONT=font, ys=75, xs=60)

status_text = widget_text(ui_outer_base, FONT=font, $
                          value = 'Ready to set up instrument.', xs=70, ys=1)

; the following base and buttons are for keyboard hotkeys
i_base=widget_base(formatbase)
i_menu=widget_button(i_base, value='', /menu)
i_open=widget_button(i_menu, value='open', resource_name='efs_open')
i_save=widget_button(i_menu, value='save', resource_name='efs_save')
i_cfg_script=widget_button(i_menu, value='cfg_script', resource_name= $
	'efs_cfg_script')
i_quit=widget_button(i_menu, value='quit', resource_name='efs_quit')
i_nod=widget_button(i_menu, value='nod', resource_name='efs_nod')
i_clear=widget_button(i_menu, value='clear', resource_name='efs_clear')
i_user=widget_button(i_menu, value='user', resource_name='efs_user')
i_oh=widget_button(i_menu, value='oh', resource_name='efs_oh')
i_neon=widget_button(i_menu, value='neon', resource_name='efs_neon')
i_argon=widget_button(i_menu, value='argon', resource_name='efs_argon')
i_krypton=widget_button(i_menu, value='krypton', resource_name='efs_krypton')
i_xenon=widget_button(i_menu, value='xenon', resource_name='efs_xenon')
i_go=widget_button(i_menu, value='go', resource_name='efs_go')
i_abort=widget_button(i_menu, value='abort', resource_name='efs_abort')
i_up=widget_button(i_menu, value='up', resource_name='efs_up')
i_down=widget_button(i_menu, value='down', resource_name='efs_down')
i_left=widget_button(i_menu, value='left', resource_name='efs_left')
i_right=widget_button(i_menu, value='right', resource_name='efs_right')
i_boxbase=widget_button(i_menu, value='', /menu)
i_newbox=widget_button(i_boxbase, value='newbox', resource_name='efs_newbox')
i_nextbox=widget_button(i_boxbase, value='nextbox', resource_name='efs_nextbox')
i_prebox=widget_button(i_boxbase, value='prebox', resource_name='efs_prebox')
i_delbox=widget_button(i_boxbase, value='delbox', resource_name='efs_delbox')


widget_control, formatbase, /realize

; initial disabled buttons
widget_control, nods, sensitive=0
widget_control, reps, sensitive=0
widget_control, obj_name, sensitive=0
widget_control, star_name, sensitive=0
widget_control, itime, sensitive=0
widget_control, stime, sensitive=0
widget_control, coadds, sensitive=0
widget_control, scoadds, sensitive=0
;widget_control, boxnum, sensitive=0
widget_control, box_button_base, sensitive=0
widget_control, i_boxbase, sensitive=0
widget_control, draw, get_value=draw_index

if serv eq  'none' then begin
    widget_control, formatmenu+7, sensitive=0
    widget_control, go_button, sensitive=0
    widget_control, abort_button, sensitive=0
endif

; unmap base that contains keyboard shortcuts
widget_control, i_base, map=0

; copy index of main draw window to a variable to be shared with plot filter
efs_draw_index=draw_index

choose_standard_filter
setup_ech_filter
; Center current box in echellogram
boxes(curbox).echphi = roundphi(echtheta)
boxes(curbox).graphi = roundphi(graphi(499))
; find where to draw box, and draw it and echellogram
fill_curbox_arrays
calc_pixel, box_x(curbox), box_y(curbox)
find_corners
draw_echellogram
draw_boxes
; tvcrs, 175, 175

xmanager, "formatbase", formatbase, /just_reg	
xmanager, 'draw', draw, /just_reg
xmanager, 'formatmenu', formatmenu, /just_reg
xmanager, 'filtermenu', filtermenu, /just_reg
; xmanager, 'thinmenu', thinmenu, /just_reg
xmanager, 'plot_filter_button', plot_filter_button, /just_reg
xmanager, 'slitmenu', slitmenu, /just_reg
xmanager, 'new_box', new_box, /just_reg
xmanager, 'del_box', del_box, /just_reg
xmanager, 'go_button', go_button, /just_reg
xmanager, 'abort_button', abort_button, /just_reg
xmanager, 'next_box', next_box, /just_reg
xmanager, 'pre_box', pre_box, /just_reg
xmanager, 'curbox_button', curbox_button, /just_reg
xmanager, 'echbox', echbox, /just_reg
xmanager, 'grabox', grabox, /just_reg
xmanager, 'itime', itime, /just_reg
xmanager, 'stime', stime, /just_reg
xmanager, 'full', full, /just_reg
xmanager, 'coadds', coadds, /just_reg
xmanager, 'scoadds', scoadds, /just_reg
xmanager, 'obj_name', obj_name, /just_reg
xmanager, 'star_name', star_name, /just_reg
xmanager, 'nods', nods, /just_reg
xmanager, 'reps', reps, /just_reg
;xmanager, 'boxnum', boxnum, /just_reg
xmanager, 'specmode', specmode, /just_reg
xmanager, 'keyboard', i_base, /just_reg
; Start xmanager to manage events
xmanager

END
;**********************************************************
