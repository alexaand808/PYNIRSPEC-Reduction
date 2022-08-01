pro browse_button_event, event
	common shared, filelist, new_line, filelist_box, add_line_box, $
		line_list, data_list, save_list_box, current_dir, save_file, $ 
		save_list_name, color_base, rgb, r, g, b, color_id, draw_id, $
		browse_button, old_list, result, displayed_list, comments, $
		add_comment_box, old_comments
	
	compare_list, old_list, list, result
	compare_list, old_comments, comments, comments_result
	result=result and comments_result
	if result eq 0 then begin
		message=['List contents have changed.', $
			'Do you wish to save list first?']
		answer=dialog_message(message, dialog_parent=event.id, $
			/question)
		if answer eq 'Yes' then begin
			save_as_button_event, event
		endif
	endif
	filelist=dialog_pickfile(file=current_dir+'/'+'user.lst', $ 
		filter='*.lst', /read, group=event.id)
	if filelist ne '' then begin
		open_list
	endif
end

pro add_line_button_event, event
	common shared

	widget_control, add_line_box, get_value=new_line
	widget_control, add_comment_box, get_value=new_comment
	if new_line[0] ne '' then begin
		if list[0] ne '' then begin
			list_size=size(list)
			if list_size(0) eq 0 then begin
				list=[list, new_line]
				comments=[comments, new_comment]
				displayed_list=string(list)+'     '+comments
			endif else begin
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
		endif else begin
			list=new_line
			comments=new_comments
			displayed_list=string(list)+'     '+comments
		endelse
	endif
	widget_control, line_list, set_value=displayed_list
end

pro remove_buttons_event, event
	common shared

	case event.value of
		0: begin
			info=widget_info(line_list, /list_select)
			;print, info, info[0]
			info_size=size(info)
			;print, 'Info size=', info_size
			temp_array=list
			temp_comments=comments
			list_size=size(list)
			;print, list_size
			if info[0] ge 0 then begin
				if info_size[0] eq 0 then begin
					j=0
					for i=0, list_size(1)-1 do begin
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
		1: begin
			warning=["This will erase all","lines from your list!"]
			answer=dialog_message(warning, /default_cancel, $
				dialog_parent=event.id, /cancel)
			if answer eq 'OK' then begin
				list=''
				comments=''
				widget_control, line_list, set_value=list
			endif
		end
		2: begin
			color_base=widget_base(/col, title='Select Color')
			color_sel=widget_draw(color_base, xs=275, ys=30)
			rgb=cw_rgbslider(color_base)
			ok_button=widget_button(color_base, value='OK')
			widget_control, color_base, /realize
			widget_control, color_sel, get_value=color_id

			xmanager, 'ok_button', ok_button, /just_reg
			xmanager, 'rgb', rgb, /just_reg
			xmanager
		end
	endcase
end

pro rgb_event, event
	common shared
	
	base = widget_info(event.id, /parent)
	if event.id eq base+1 then r=event.value
	if event.id eq base+2 then g=event.value
	if event.id eq base+3 then b=event.value
	wset, color_id
	color=r+256L*(g+256L*b)
	polyfill, [0,1,1,0], [0,0,1,1], $
		color=color
end

pro ok_button_event, event
	common shared

	wset, draw_id
	color=r+256L*(g+256L*b)
	polyfill, [0,1,1,0], [0,0,1,1], color=color	
	widget_control, color_base, /destroy
end
	
pro save_as_button_event, event
	common shared

	widget_control, save_list_box, get_value=save_list_name
	if strmid(save_list_name[0], 0, 1) eq '/' then begin
		default=save_list_name[0]
	endif else begin
		default=current_dir+'/'+save_list_name[0]
	endelse
	save_file=dialog_pickfile(/write,  file=default, filter='*.lst', $
		group=event.id)
	print, save_file
	if save_file ne '' then begin
		result=findfile(save_file)
		if result[0] eq '' then begin
			save_list
		endif else begin
			message=['File ' + save_file+' exists.', $
				'Do you wish to Overwrite?']
			answer=dialog_message(message,title='Overwrite File?',$
				/default_no, dialog_parent=event.id, /question)
			if answer eq 'Yes' then save_list
		endelse
	endif
end

pro control_buttons_event, event
	common shared

	case event.value of
		0: begin
		compare_list, old_list, list, result
		compare_list, old_comments, comments, comments_result
		result=result and comments_result
		if result eq 0 then begin
			message=['List contents have changed.', $
				'Do you wish to save list first?']
			answer=dialog_message(message, dialog_parent=event.id,$
				/question)
			if answer eq 'Yes' then begin
				save_as_button_event, event
			endif
		endif
		widget_control, event.top, /destroy
		end
		2: widget_control, event.top, /destroy
		else: print, event.value
	endcase
end

pro open_list
	common shared

	flt_array=fltarr(1000)
	temp_comments=strarr(1000)
	read_lines, filelist, flt_array, temp_comments, i, err 
	if err eq -215 then begin
		message=['File does not exist', $
			 'Do you wish to create a new one?']
		answer=dialog_message(message, dialog_parent=browse_button, $
			/question)
		if answer eq 'Yes' then begin
			list=''
			widget_control, filelist_box, set_value=filelist
			widget_control, save_list_box, set_value=filelist
		endif
	endif

	if i ge 1 then begin
		list=fltarr(i)
		comments=strarr(i)
		for j=0, i-1 do begin
			list(j)=flt_array(j)
			comments(j)=temp_comments(j)
		endfor
	endif else begin
		list=''
		comments=''
	endelse

	widget_control, filelist_box, set_value=filelist
	widget_control, save_list_box, set_value=filelist

	old_list=list
	old_comments=comments
	displayed_list=string(list)+'     '+comments
	widget_control, line_list, set_value=displayed_list
end

pro save_list
	common shared

	list_size=size(list)
	openw, 1, save_file, error=err
	print, err
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
	endif else begin
		if list[0] eq '' then begin
			printf, 1, ''
		endif else begin
			if list_size[0] eq 0 then begin
				printf, 1, list+' '+comments
			endif else begin
				for i=0, list_size(1)-1 do begin
					line=string(list(i))+' '+comments(i)
					printf, 1, line
				endfor
			endelse
		endelse
		close, 1
		old_list=list
		old_comments=comments
	endelse
end

pro user_lines, list, color
	common shared

print, 'in user_lines.pro... ', ' list=',list, ' color=', color
filelist='Click browse to open file...'
; list=[2.0,4.0,5.0,8.0,9.0,34.0,56.0]
comments=''
old_list=list
old_comments=comments
new_line=''
remove=['Remove Selection', 'Clear List', 'Select Color']
control=['OK', 'Apply', 'Cancel']
cd, '.', current=current_dir
r=0 & g=0 & b=0
result=0

base=widget_base(/col, Title='User Line List')
file_base=widget_base(base, /row)
filelist_box=cw_field(file_base, value=filelist, title='List Filename:', $
	xs=30, /noedit)
browse_button=widget_button(file_base, value='Browse')
line_list=widget_list(base, value=string(list), ys=4)
add_line_base=widget_base(base, /row)
add_line_box=cw_field(add_line_base, value=new_line, /floating, $
	title='Add Line to List:', xs=30)
add_line_button=widget_button(add_line_base, value='Add')
add_comment_box=cw_field(base, value='', title='New Line Comment:', xs=35)
remove_base=widget_base(base, /row)
remove_buttons=cw_bgroup(remove_base, remove, /row)
draw_base=widget_base(remove_base, /col)
color_draw=widget_draw(draw_base, xs=60, ys=20)
save_base=widget_base(base, /row)
save_list_box=cw_field(save_base, value='No File Opened...', xs=30, $
	title='Save List As:') 
save_as_button=widget_button(save_base, value='Save')
control_buttons=cw_bgroup(base, control, /row)

widget_control, base, /realize
widget_control, color_draw, get_value=draw_id

xmanager, 'browse_button', browse_button, /just_reg, /no_block
xmanager, 'add_line_button', add_line_button, /just_reg, /no_block
xmanager, 'remove_buttons', remove_buttons, /just_reg, /no_block
xmanager, 'save_as_button', save_as_button, /just_reg, /no_block
xmanager, 'control_buttons', control_buttons, /just_reg, /no_block
xmanager

end


































































































