;+
; NAME:
;	CW_TOGGLE
;
; PURPOSE:
;	This widget consists of a number of buttons, of which, at most
;	one can be selected at a time.  The widget consists of a label
;	and a draw widget in which buttons are drawn in.  CW_TOGGLE's
;	value is an integer ranging from 0 to the number of buttons
;	minus one, representing the index of the button pressed
;	
;	
; CATEGORY:
;	Widget Clusters.
;
; CALLING SEQUENCE:
;	Result = CW_TOGGLE(Parent, Names)
;
; INPUTS:
;	Parent:	The widget ID of the widget to be the field's parent.
;
;        Names: An array of strings to give the labels of the buttons.
;               The number of buttons is also determined by this
;               array.
;
; KEYWORD PARAMETERS:
;	TITLE:	A string containing the text to be used as the label for the
;		field.  If omitted, not label is written.
;
;	VALUE:	The initial value of the widget.  This value is
;		the index of the button to be pressed.  It is zero
;		based, starting with 0 as the left-most button.  The
;		default is 0.
;
;      UVALUE:	A user value to assign to the draw cluster.  This value
;		can be of any type.
;
;	XSIZE:	An explicit horizontal size (in pixels) for the draw widget
;		area.  The default size is 120.
;
;	YSIZE:	An explicit vertical size (in pixels) for the draw widget
;		area.  The default is 30.
;
;	 FONT:	A string containing the name of the X Windows font to use
;		for the TITLE of the field.
;
;   EVENT_PRO:   A string containing the name of the event handler for
;               the toggle widget event.
;
; OUTPUTS:
;	This function returns the widget ID of the newly-created cluster.
;
; COMMON BLOCKS:
;	None.
;
; PROCEDURE:
;	Create the widgets, set up the appropriate event handlers, and return
;	the widget ID of the newly-created cluster.
;
; EXAMPLE:
;	The code below creates a main base with a toggle cluster attached
;	to it.  The cluster has the title "LIGHT:", two buttons
;	labeled 'On' and 'Off', with 'Off' set as the initial setting.
;
;		base = WIDGET_BASE()
;		field = CW_TOGGLE(base, ['On', 'Off'], TITLE="LIGHT:", value=1)
;		WIDGET_CONTROL, base, /REALIZE
;
; MODIFICATION HISTORY:
; 	Written by: Jason L. Weiss, January 1999
;    February 1999: JLW -- Added support for more than two buttons
;-

function col, r, g, b
; convert r, g, b values to a color
; syntax, color=col(r,g,b)

	color=r+(256L*(g+(256L*b)))
	return, color

end

function cw_toggle_event, event
; main handler for hitting the toggle

; if button was clicked
if event.type eq 0 then begin
        
        ; get uval
	base=event.handler
	stash=widget_info(base, /child)
	widget_control, stash, get_uval=state
	ys=state.ys
	xs=state.xs
	num=state.num	
	fnum=float(num)

        ; new_val is the index of the part of the toggle clicked
	new_val=event.x/(xs/num)
        ; if it is a new value, then redraw buttons
	if state.value ne new_val then begin

            save = !d.window
            wset, state.draw_index
            state.value=new_val

            ; draw buttons
            for i=0, num-1 do begin

                ; set colors
                if i eq state.value then begin
                    base_color=state.gray.norm
                    top_edge=state.shade.dark
                    bot_edge=state.shade.light
                endif else begin
                    base_color=state.gray.light
                    top_edge=state.shade.light
                    bot_edge=state.shade.dark
                endelse

                ; draw buttons in
                polyfill, [i/fnum, (i+1)/fnum, (i+1)/fnum, i/fnum], $
                  [0, 0, 1, 1], color=base_color, /line_fill, /normal

                ; draw edges of buttons for 3D look
                plots, (i/fnum)*xs, 0, /device
                plots, (i/fnum)*xs, ys-1, color=top_edge, /continue, thick=1, $
                  /device
                plots, (xs*(i+1)/fnum)-1, ys-1, color=top_edge, /continue, $
                  thick=1, /device
                plots, (xs*(i+1)/fnum)-1, 0, color=bot_edge, /continue, $
                  thick=1, /device
                plots, (i/fnum)*xs, 0, color=bot_edge, /continue, thick=1, $
                  /device
                plots, (xs*(i/fnum))+1, 1, /device
                plots, (xs*(i/fnum))+1, ys-2, color=top_edge, /continue, $
                  thick=1, /device
                plots, (xs*(i+1)/fnum)-2, ys-2, color=top_edge, /continue, $
                  thick=1, /device
                plots, (xs*(i+1)/fnum)-2, 1, color=bot_edge, /continue, $
                  thick=1, /device
                plots, (xs*(i/fnum))+2, 1, color=bot_edge, /continue, $
                  thick=1, /device
                
                ; write out names
                xyouts, (xs*(i+0.5)/fnum), (ys/2)-5, /device, state.names(i), $
                  color=0, font=0, charthick=2, charsize=1.5, alignment=0.5

            endfor

            wset, save
                ; send event to event handler
		epro=state.epro	
		ret = {ID:state.base, TOP:event.top, HANDLER:0L, $
			value:state.value}

		widget_control, stash, set_uval=state, /no_copy
		
		if epro ne '' then CALL_PROCEDURE, epro, ret
		ret=event
	endif else begin
		ret=0
	endelse
endif else begin
	ret=0
endelse

return, ret

end

function cw_toggle_get_value, id
; handles widget_control, get_value call

        ; get value from uval and return it
	stash=widget_info(id, /child)
	widget_control, stash, get_uvalue=state, /no_copy
	ret=state.value
	widget_control, stash, set_uvalue=state
	return, ret
end

pro cw_toggle_set_value, id, value
; handles widget_control, set_value call

	; get uval
	stash=widget_info(id, /child)
	widget_control, stash, get_uvalue=state, /no_copy
        wset, state.draw_index
        ; set state.value to a value other than desired so that event will
        ; be passed and applied
	state.value= not(value)
        ; set clicked point as middle of button
	x=state.xs*(value+0.5)/float(state.num)
        ; define fake event
	event={ID:id, top:id, handler:id, type:0, x:x, y:state.ys/2, $
               press:1, release:0, clicks:1}
        ; send it to event handler
	widget_control, id, send_event=event
	widget_control, stash, set_uvalue=state

end

pro cw_toggle_realize, id
; called when toggle is created.  basically just fills the draw widget
; with the buttons

; get uval
stash=widget_info(id, /child)
widget_control, stash, get_uval=state, /no_copy
; get draw widget index
widget_control, state.draw_id, get_value=draw_index
state.draw_index=draw_index
widget_control, stash, set_uval=state
; get draw size
ys=state.ys
xs=state.xs
; get number of buttons
num=state.num
fnum=float(num)
save=!d.window
wset, draw_index

; draw buttons
for i=0, num-1 do begin

; assign colors
if i eq state.value then begin
	base_color=state.gray.norm
	top_edge=state.shade.dark
	bot_edge=state.shade.light
endif else begin
	base_color=state.gray.light
	top_edge=state.shade.light
	bot_edge=state.shade.dark
endelse

; fill buttons
polyfill, [i/fnum, (i+1)/fnum, (i+1)/fnum, i/fnum], [0, 0, 1, 1], $
	color=base_color

; draw edges of buttons for 3D look
plots, (i/fnum)*xs, 0, /device
plots, (i/fnum)*xs, ys-1, color=top_edge, /continue, thick=1, /device
plots, (xs*(i+1)/fnum)-1, ys-1, color=top_edge, /continue, thick=1, /device
plots, (xs*(i+1)/fnum)-1, 0, color=bot_edge, /continue, thick=1, /device
plots, (i/fnum)*xs, 0, color=bot_edge, /continue, thick=1, /device
plots, (xs*(i/fnum))+1, 1, /device
plots, (xs*(i/fnum))+1, ys-2, color=top_edge, /continue, thick=1, /device
plots, (xs*(i+1)/fnum)-2, ys-2, color=top_edge, /continue, thick=1, /device
plots, (xs*(i+1)/fnum)-2, 1, color=bot_edge, /continue, thick=1, /device
plots, (xs*(i/fnum))+2, 1, color=bot_edge, /continue, thick=1, /device

; fill buttons with text
xyouts, (xs*(i+0.5)/fnum), (ys/2)-5, /device, state.names(i), color=0, $
	font=0, charthick=2, charsize=1.5, alignment=0.5

endfor
wset, save

end

function cw_toggle, parent, names, title=title, value=value, xs=xs, ys=ys, $
	uvalue=uvalue, font=font, event_pro=epro
; main program that sets up widgets and uvals

; define colors
red={light:col(255, 127, 127), norm:255, dark:200}
gray={light:col(192, 192, 192), norm:col(170, 170, 170), dark:col(127, 127, 127)}
shade={light:col(235, 235, 235), dark:col(75, 75, 75)}
num=n_elements(names)

if not keyword_set(value) then value=0
if not keyword_set(xs) then xs=120
if not keyword_set(ys) then ys=30
if not keyword_set(uvalue) then uvalue=0
if not keyword_set(epro) then epro=''

xs=xs>60
ys=ys>20
on_error, 2

state=({value:value, base:0, draw_id:0, draw_index:0, red:red, gray:gray, $
	xs:xs, ys:ys, shade:shade, names:names, epro:epro, num:num})

base=widget_base(parent, event_func='cw_toggle_event', uvalue=uvalue, $
	func_get_value='cw_toggle_get_value', $
	pro_set_value='cw_toggle_set_value', $
	notify_realize='cw_toggle_realize')
if keyword_set(title) then begin
	label_base=widget_base(base, /row)
    if keyword_set(font) then begin
	label=widget_label(label_base, font=font, value=title)
    endif else begin
	label=widget_label(label_base, value=title)
    endelse
	draw=widget_draw(label_base, xs=xs, ys=ys, /button_events)
endif else begin
	draw=widget_draw(base, xs=xs, ys=ys, /button_events)
endelse

state.base=base
state.draw_id=draw
widget_control, widget_info(base, /child), set_uvalue=state, /no_copy

return, base

end
