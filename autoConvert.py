import win32api
import win32con
import os

for i in range(23,61):
#	filename = '(' + str(i) + ')'

	filename = 'flake' + str(i) + '.' + 'exr'
	os.system(filename)

	# win32api.keybd_event(17,0,0,0)
	# win32api.keybd_event(83,0,0,0)
	# win32api.keybd_event(83,0,win32con.KEYEVENTF_KEYUP,0)
	# win32api.keybd_event(17,0,win32con.KEYEVENTF_KEYUP,0)
	# num = int(i/10)+48
	# num2 = i%10+48
	# win32api.keybd_event(num,0,0,0)
	# win32api.keybd_event(num2,0,0,0)
	# win32api.keybd_event(num,0,win32con.KEYEVENTF_KEYUP,0)
	# win32api.keybd_event(num2,0,win32con.KEYEVENTF_KEYUP,0)
	# win32api.keybd_event(13,0,0,0)
	# win32api.keybd_event(13,0,win32con.KEYEVENTF_KEYUP,0)
	# win32api.keybd_event(17,0,0,0)
	# win32api.keybd_event(87,0,0,0)
	# win32api.keybd_event(87,0,win32con.KEYEVENTF_KEYUP,0)
	# win32api.keybd_event(17,0,win32con.KEYEVENTF_KEYUP,0)


