; plot filter.pro 
;
; this procedure plots a profile of the throughput of the filters used
; in nirspec
;
;  written by Jason L. Weiss

pro quit_event, event
        common shared, plot_draw, base, type, num, draw_index, ech, gra, $
          data_x, data_y, data_ech, data_gra, graorder, quit, efs_draw_id, $
          atm_x, atm_y, index, plot_atm_y, temp_atm_x, temp_atm_y, $
          print_flag, printer_name, font, printname_box

        common plot_shared

        ; set plotted flag to zero since window is destroyed
        plotted=0
        widget_control, base, /destroy

end

pro print_buttons_event, event
        common shared


        case event.value of
            'OK': begin
                widget_control, printname_box, get_value=temp
                printer_name=temp[0]
                draw_profile
                widget_control, event.top, /destroy
            end
            'CANCEL': begin
                print_flag=0
                widget_control, event.top, /destroy
            end
        endcase
end

pro print_event, event
        common shared


        if print_flag eq 0 then begin
            print_flag=1
            print_base=widget_base(/col, title='Print')
            printname_box=cw_field(print_base, value=printer_name, $
                                   title='Printer Name: ', font=font)
            print_buttons=cw_bgroup(print_base, ['OK', 'CANCEL'], font=font, $
                                   /return_name, /row)
            
            widget_control, print_base, /realize, group_leader=base
            xmanager, 'print_buttons', print_buttons, /just_reg, /no_block
            xmanager
        endif

end


pro base_event, event
        common shared
        common plot_shared

        ; resize events
        xs=event.x
        ys=event.y

        widget_control, plot_draw, xs=xs
        widget_control, plot_draw, ys=ys

        ; redraw profile
        draw_profile

end

pro base_killed, event
        common shared
        common plot_shared
        
        ; set plotted flag to zero to note that window was destroyed
        plotted=0
        wset, efs_draw_id
end

pro draw_profile
        common shared
        common shared_3
        common plot_shared

        red=255
        white=255+256L*(255+(255*256L))
        green=0+256L*(255+(0*256L))
        blue=0+256L*(0+(255*256L))
	black=0

        wset, draw_index

        ; calculate the blaze angle
        echblazelam=ech.sigma*cos(ech.gamma)*2.0*sin(ech.theta)
;        gralambdanot=2.0*gra.sigma*sin(gra.theta)/graorder
;        gralambdanot=gra.sigma*(sin(gra.haang+graphi(499))+ $
;                                sin(graphi(499)-gra.haang))/graorder
        gralambdanot=gra.sigma*(sin(0.61086)+ $
                                sin(-0.2618))/graorder

        ;for which filter are we displaying the plot
        case type of
            'NIRSPEC-1': begin
                data_file='data/nir1data.txt'
                atm_file='data/nir1atm.dat'
                end
            'NIRSPEC-2': begin
                data_file='data/nir2data.txt'
                atm_file='data/nir2atm.dat'
                end
            'NIRSPEC-3(J)': begin
                data_file='data/nir3data.txt'
                atm_file='data/nir3atm.dat'
                end
            'NIRSPEC-4': begin
                data_file='data/nir4data.txt'
                atm_file='data/nir4atm.dat'
                end
            'NIRSPEC-5(3rd)': begin
                data_file='data/nir5data.txt'
                atm_file='data/nir5atm.dat'
                end
            'NIRSPEC-6': begin
                data_file='data/nir6data.txt'
                atm_file='data/nir6atm.dat'
                end
            'NIRSPEC-7': begin
                data_file='data/nir7data.txt'
                atm_file='data/nir7atm.dat'
                end
            'BR-GAMMA(2.16)': begin
                data_file='data/bgamdata.txt'
                atm_file='data/bgamatm.dat'
                end
            'CO(2.3)': begin
                data_file='data/codata.txt'
                atm_file='data/coatm.dat'
                end
            'FeII(1.64)': begin
                data_file='data/feiidata.txt'
                atm_file='data/feiiatm.dat'
                end
            'H2(2.12)': begin
                data_file='data/h2data.txt'
                atm_file='data/h2atm.dat'
                end
            'K': begin
                data_file='data/kdata.txt'
                atm_file='data/katm.dat'
                end
            'K-PRIME': begin
                data_file='data/kprimedata.txt'
                atm_file='data/kprimeatm.dat'
                end
            'L-PRIME': begin
                data_file='data/lprimedata.txt'
                atm_file='data/lprimeatm.dat'
                end
            'M-PRIME': begin
                data_file='data/mprimedata.txt'
                atm_file='data/mprimeatm.dat'
                end
            'KL': begin
                data_file='data/kldata.txt'
                atm_file='data/klatm.dat'
                end
            'M-WIDE': begin
                data_file='data/mwidedata.txt'
                atm_file='data/mwideatm.dat'
                end
            'Thin Blocker': begin
                data_file='data/pk50thindata.txt'
                atm_file='data/pk50thinatm.dat'
                end
            'Thick Blocker': begin
                data_file='data/pk50thckdata.txt'
                atm_file='data/pk50thckatm.dat'
                end
            else: print, "error"
        endcase
        
                ; get data from the file
                read_data, data_file, data_x, data_y, num

                ; get atmospheric transmission curve
                min=min(data_x[0:num-1], max=max)
                min=min/1000.0
                max=max/1000.0
                print, min, max
                atmnum=0
                read_data, atm_file, temp_atm_x, temp_atm_y, atmnum
                temp_x=data_x/1000.0

                plot_atm_y=interpol(temp_atm_y[0:atmnum-1], $
                                    temp_atm_x[0:atmnum-1], temp_x[0:num-1])
                
                ; calculate and overplot the reduction in throughput
                ; due to the echell grating and the cross disperser
                for i=0, num-1 do begin 
                    ;index=where((atm_x ge ((data_x[i]-0.2)/1000.0)) and $
                    ;            (atm_x le ((data_x[i]+0.2)/1000.0)))
                    ;if index[0] gt 0 then $
                    ;  plot_atm_y[i]=float(atm_y[index[0]]) $
                    ;else plot_atm_y[i]=1
                    data_y[i]=float(data_y[i])*float(plot_atm_y[i])
                    echorder=fix(echblazelam*1000000.0/data_x(i)+0.5)
                    gragamma=(graorder*!pi*1000000.0*(gralambdanot- $
                             data_x(i)/1000000.0))/data_x(i)

;                    gragamma=(!pi*100000.0*gra.sigma*(sin(gra.haang+$
;                    graphi(499))-sin(gra.haang-graphi(499))))/(data_x[i])
;                    print, gragamma
                    data_gra(i)=(((sin(gragamma))/gragamma)^2)*data_y(i)
                    echgamma=(echorder*!pi*1000000.0*(echblazelam/echorder- $
                          data_x(i)/1000000.0))/data_x(i)
;                    echgamma=(echorder*!pi*1000000.0*(echblazelam/(echorder* $
;                                                      data_x(i))))
                    data_ech(i)=(((sin(echgamma))/echgamma)^2)*data_gra(i)
                    plot_atm_y[i]=100.0*plot_atm_y[i]
                endfor


                if print_flag eq 1 then begin
                    set_plot, 'ps'
                    device, filename='idl.ps'
		    black=0 & red=100 & green=150 & blue=200 & white=255
                endif
                
                plot, data_x[0:num-1], plot_atm_y, color=black, $
                  title=type, ymargin=[7,2], xtitle='Wavelength (nm)', $
                  ytitle='% Transmission', yrange=[0, 100], background=white,/nodata
		oplot, data_x[0:num-1], plot_atm_y, color=blue
                oplot, data_x(0:num-1), data_ech, color=red
                oplot, data_x(0:num-1), data_gra, color=black
                oplot, data_x(0:num-1), data_y, color=green

                ; plot legend
                window_size=widget_info(plot_draw, /geometry)
                left=60
                right=window_size.xsize-18
                ymid=16
                lmid=left+((right-left)/3.0)-24
                rmid=left+((right-left)*2.0/3.0)-88
		if print_flag eq 1 then goto,skipout
                plots, [left, left+24], [ymid+3, ymid+3], /device, color=blue
                plots, [lmid-30, lmid-6], [ymid+3, ymid+3], /device, color=green
                plots, [rmid-30, rmid-6], [ymid+3, ymid+3], /device, color=black
                plots, [right-128, right-104], [ymid+3, ymid+3], /device, $
                  color=red
                xyouts, left+30, ymid, 'Atmosphere', /device, alignment=0.0, color=blue
                xyouts, lmid, ymid, 'With Filter', /device, color=green, $
                  alignment=0.0
                xyouts, rmid, ymid, 'With Cross Dispersor', /device, color=black, $
                  alignment=0.0
                xyouts, right, ymid, 'With Echelle Grating', /device, color=red, $
                  alignment=1.0
		
		skipout:

                if print_flag eq 1 then begin
                    device, /close
		    set_plot,'x'
                    command=strcompress('lpr -P'+printer_name+' idl.ps')
                    spawn, command, result
                    if result(0) ne '' then begin
                        answer=dialog_message(result, dialog_parent=base, /error)
                    endif
                    spawn, '\rm idl.ps'
                    print_flag=0
                endif

                wset, efs_draw_id
                draw_echellogram
                draw_boxes
end

pro plot_filter, filter, echinfo, grainfo
        common shared
        common plot_shared

        ; rename local variables
        ech=echinfo
        gra=grainfo
        print, ech
        print, gra
        efs_draw_id=efs_draw_index
	font='-adobe-times-bold-r-normal--20-140-100-100-p-100-iso8859-1'

        type=filter.name
        graorder=filter.graorder
        num=0
        print_flag=0
        printer_name=''
        ; filter arrays
        maxlines=5000
        atmmax=40000
        data_x=fltarr(maxlines)
        data_y=fltarr(maxlines)
        index=fltarr(atmmax)
        atm_x=fltarr(atmmax)
        atm_y=fltarr(atmmax)
        temp_atm_x=fltarr(atmmax)
        temp_atm_y=fltarr(atmmax)
        plot_atm_y=fltarr(maxlines)
        data_ech=fltarr(maxlines)
        data_gra=fltarr(maxlines)
        data_atm=fltarr(maxlines)

        ; draw window and widgets if not already drawn
        if plotted eq 0 then begin
            ;set up widgets
            base=widget_base(/col, title="Filter Profile", /TLB_SIZE_EVENTS, $
                            kill_notify='base_killed', group_leader=efs_base)
            plot_draw=widget_draw(base, xs=640, ys=480)
            print=widget_button(base, value='Print', font=font)
            quit=widget_button(base, value='Quit', font=font)
        
            widget_control, base, /realize
            widget_control, plot_draw, get_value=draw_index        
        endif

        ; set plotted flag to one since window is drawn
        plotted=1

	draw_profile
        
        xmanager, 'base', base, /just_reg, /no_block
        xmanager, 'quit', quit, /just_reg, /no_block
        xmanager, 'print', print, /just_reg, /no_block
        xmanager

        print, 'done'
end
