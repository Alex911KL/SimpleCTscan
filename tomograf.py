import math
import skimage as sk
import numpy as np
import datetime
import re
import tkinter as tk
from tkinter import ttk, filedialog
from ttkthemes import ThemedTk
import pydicom as dic
from datetime import datetime
from PIL import Image, ImageTk


# Normalizacja obrazu do zakresu [0, 255]
def rescale(img):
  normalized = np.zeros(shape=(img.shape))
  img_max = img.max()
  img_min = 0

  for x in range(img.shape[0]):
    normalized[x] = np.interp(img[x], (img_min, img_max), (0, 255))
  return normalized


# Załadowanie obrazu
def loadImage(path):
  img = sk.io.imread(path, as_gray=True)
  img = rescale(img)
  return img


# Algorytm Brasenhama - zwraca średnią jasność
def brasenhamAlgorithm(img, x1, y1, x2, y2):
  max_x, max_y = len(img), len(img[0])
  sum_brightness, n = 0, 0
  dx, dy = abs(x2 - x1), abs(y2 - y1)
  x_inc, y_inc = (1 if x1 < x2 else -1), (1 if y1 < y2 else -1)
  x, y = x1, y1

# Oś wiodąca OX (większa szybkość zmian)
  if dx > dy:
    err, d_err = 2 * dy - dx, 2 * dy
    for _ in range(dx):
      if 0 <= x < max_x and 0 <= y < max_y:
        sum_brightness += img[x][y]
        n += 1
      if err >= 0:
        y += y_inc
        err -= 2 * dx
      err += d_err
      x += x_inc
  else:
    err, d_err = 2 * dx - dy, 2 * dx
    for _ in range(dy):
      if 0 <= x < max_x and 0 <= y < max_y:
        sum_brightness += img[x][y]
        n += 1
      if err >= 0:
        x += x_inc
        err -= 2 * dy
      err += d_err
      y += y_inc

  return sum_brightness / n if n else 0

# Odwrotny algorytm Brasenhama - zwraca koordynaty
def inverseBresenhamAlgorithm(max_x, max_y, x1, y1, x2, y2):
  coords = []
  dx, dy = abs(x2 - x1), abs(y2 - y1)
  x_inc, y_inc = (1 if x1 < x2 else -1), (1 if y1 < y2 else -1)
  x, y = x1, y1

  if dx > dy:
    err, d_err = 2 * dy - dx, 2 * dy
    for _ in range(dx):
      if 0 <= x < max_x and 0 <= y < max_y:
        coords.append((x, y))
      if err >= 0:
        y += y_inc
        err -= 2 * dx
      err += d_err
      x += x_inc
  else:
    err, d_err = 2 * dx - dy, 2 * dx
    for _ in range(dy):
      if 0 <= x < max_x and 0 <= y < max_y:
        coords.append((x, y))
      if err >= 0:
        x += x_inc
        err -= 2 * dy
      err += d_err
      y += y_inc

  return coords

# Transformata Radona
def radonTransform(img, img_label, emitterRange = 180, numOfDetectors = 180, numOfScans = 180, step = 2):
  label_progress("Transformata Radona")
  step = math.radians(step)
  center = (len(img)//2, len(img[0])//2)
  R = max(center) * math.sqrt(2)

  alpha = math.pi/2
  phi = math.pi*emitterRange/180
  sinogram = np.zeros(shape=(numOfScans,numOfDetectors))

# Pętla do skanów
  for scan in range(numOfScans): # dla każdego skanu
    progress((scan / numOfScans) * 100)
    xe = center[0] + round(R * math.cos(alpha))
    ye = center[1] - round(R * math.sin(alpha))

    for det in range(numOfDetectors):
      xd = center[0] + round(R * math.cos(alpha+math.pi-phi/2 + det*phi/(numOfDetectors-1)))
      yd = center[1] - round(R * math.sin(alpha+math.pi-phi/2 + det*phi/(numOfDetectors-1)))
      sinogram[scan][det] = brasenhamAlgorithm(img, xe, ye, xd, yd) #Linia od emitera do detektora

    alpha += step

  showImage(rescale(sinogram), img_label)
  return rescale(sinogram)

#Filr z szarym tłem
def filtrGREY(sinogram, img_label):
  label_progress("Filtrowanie")

  sinogram_norm = sinogram / 255.0

  #Filtr Ram-Lak
  n = 10  # Połowa maski
  kernel = np.array([(-4 / (math.pi ** 2 * i ** 2)) if i % 2 else 0 for i in range(1, n + 1)])
  kernel = np.concatenate((kernel[::-1], [1], kernel))  #Pełna maska

  sinogram_filtered = np.zeros_like(sinogram_norm)
  for i in range(sinogram.shape[0]):
    sinogram_filtered[i, :] = np.convolve(sinogram_norm[i, :], kernel, mode="same")

  sinogram_filtered -= np.min(sinogram_filtered)
  if np.max(sinogram_filtered) > 0:
    sinogram_filtered = 255 * (sinogram_filtered / np.max(sinogram_filtered))

  showImage(sinogram_filtered.astype(np.uint8), img_label, flag=True)
  return sinogram_filtered

def filtr(sinogram, img_label):
  label_progress("Filtrowanie")

  #Filtr Ram-Lak
  n = 10  # Połowa maski
  kernel = np.array([(-4 / (math.pi ** 2 * i ** 2)) if i % 2 else 0 for i in range(1, n + 1)])
  kernel = np.concatenate((kernel[::-1], [1], kernel))  # Pełna maska

  sinogram_filtered = sinogram.copy()

  for i in range(sinogram.shape[0]):
    progress((i / sinogram.shape[0]) * 100)
    sinogram_filtered[i, :] = np.convolve(sinogram[i, :], kernel, mode="same")

#Ustawianie wartości na 0, gdy < 0
  if np.min(sinogram_filtered) < 0:
    for i in range(sinogram_filtered.shape[0]):
      for j in range(sinogram_filtered.shape[1]):
        if sinogram_filtered[i, j] < 0:
          sinogram_filtered[i, j] = 0

  showImage(rescale(sinogram_filtered), img_label, flag=True)
  return sinogram_filtered


# Odwrotna transformata Radona
def inverseRadonTransform(sinogram, img,img_label, emitterRange = 180, numOfDetectors = 180, numOfScans = 180, step = 2):

  label_progress("Odwrotna Transformata Radona")
  step = math.radians(step)
  center = (img.shape[0]//2, img.shape[1]//2)
  R = max(center) * math.sqrt(2)
  alpha = math.pi/2
  phi = math.pi*emitterRange/180
  reconstructed = np.zeros(shape=img.shape)

  for scan in range(numOfScans):
    progress((scan/numOfScans)*100)
    xe = center[0] + round(R * math.cos(alpha))
    ye = center[1] - round(R * math.sin(alpha))

    for det in range(numOfDetectors):
      xd = center[0] + round(R * math.cos(alpha+math.pi-phi/2 + det*phi/(numOfDetectors-1)))
      yd = center[1] - round(R * math.sin(alpha+math.pi-phi/2 + det*phi/(numOfDetectors-1)))
      coords = inverseBresenhamAlgorithm(img.shape[0], img.shape[1], xe, ye, xd, yd) #Linia między emiterem a detektorem
      reconstructed[tuple(np.transpose(coords))] += sinogram[scan][det] #Wzmacniamy pixele wyznaczonej linii o odpowiednią średnią
    alpha += step
    mse, rmse = RMSE(img, rescale(reconstructed))
    label_progress(f"Odwrotna Transformata Radona\nRMSE: {rmse:.4f}")


  showImage(rescale(reconstructed), img_label)
  return rescale(reconstructed)

#Obliczanie błędu średnio-kwadratowego
def RMSE(img1, img2):
  mse = np.square(np.subtract(img1, img2)).mean()
  return (mse, math.sqrt(mse)) # (mse, rmse)

#Zapisywanie do DICOM
def convertToDCM(reconstructed, name="BRAK", patient_id="0", date="", comment="BRAK"):

  #Meta tagi
  file_meta = dic.dataset.FileMetaDataset()
  file_meta.MediaStorageSOPClassUID = dic.uid.UID('1.2.840.10008.5.1.4.1.1.2')
  file_meta.MediaStorageSOPInstanceUID = dic.uid.generate_uid()
  file_meta.ImplementationClassUID = dic.uid.UID("1.2.826.0.1.3680043.8.498.1")
  dcm = dic.dataset.FileDataset("output.dcm", {}, file_meta=file_meta, preamble=b"\0" * 128)
  dcm.file_meta.TransferSyntaxUID = dic.uid.ImplicitVRLittleEndian
  dcm.is_little_endian = True
  dcm.is_implicit_VR = True

  #Dane o pacjencie i badaniu
  dt = datetime.now()
  if date == "":
    dcm.StudyDate = dt.strftime('%Y%m%d')
    dcm.ContentDate = dt.strftime('%Y%m%d')
  else:
    date = f"{date[6:10]}{date[3:5]}{date[0:2]}"
    dcm.StudyDate = date
    dcm.ContentDate = date
  dcm.PatientName = name
  dcm.PatientID = patient_id
  dcm.StudyID = "1234"
  dcm.SeriesNumber = "1"
  dcm.PatientComments = comment

  #Unikatowe ID
  dcm.SOPInstanceUID = dic.uid.generate_uid()
  dcm.SeriesInstanceUID = dic.uid.generate_uid()
  dcm.StudyInstanceUID = dic.uid.generate_uid()
  dcm.FrameOfReferenceUID = dic.uid.generate_uid()

  #Parametry obrazu
  dcm.ImageType = ["ORIGINAL", "PRIMARY", "AXIAL"]
  dcm.Modality = "CT"
  dcm.Rows = reconstructed.shape[0]
  dcm.Columns = reconstructed.shape[1]
  dcm.BitsAllocated = 8
  dcm.BitsStored = 8
  dcm.HighBit = dcm.BitsStored - 1
  dcm.SamplesPerPixel = 1
  dcm.PhotometricInterpretation = 'MONOCHROME2'
  dcm.PixelRepresentation = 0
  dcm.PixelData = reconstructed.astype(np.uint8).tobytes() # piksele obrazka

  #Zapisanie własnych danych
  block = dcm.private_block(0x000b, "PUT 155932 155885", create=True)
  block.add_new(0x01, "SH", comment)
  dcm.save_as("output.dcm", write_like_original=False)






############
#### UI ####
############

#
def label_progress(val):
  progress_bar_label.config(text=val)
  root.update()

def toggle_dicom_fields():
  if dicom_fields_frame.winfo_viewable():
    dicom_fields_frame.pack_forget()
  else:
    dicom_fields_frame.pack()


def progress(val):
  if progress_bar['value'] < 100:
    progress_bar['value'] = val
  root.update()

#Wybór obrazu
def choose_file():
  global file_path, im
  name_entry.delete(0, tk.END)
  id_entry.delete(0, tk.END)
  date_entry.delete(0, tk.END)
  comment_entry.delete(0, tk.END)
  im = None
  file_path = filedialog.askopenfilename(filetypes=[("Pliki", ".jpg .png .dcm")])
  file_label.config(text=file_path)

  if file_path.endswith(".dcm"):
    try:
      dcm = dic.dcmread(file_path)
      im = rescale(dcm.pixel_array)

      #Automatyczne zaznaczenie checkboxa "Wygeneruj DICOM"
      dicom_var.set(1)
      dicom_fields_frame.pack()

      #Wypełnienie pól danymi DICOM
      name_entry.insert(0, str(dcm.get('PatientName', 'Brak danych')))
      id_entry.insert(0, str(dcm.get('PatientID', 'Brak danych')))
      study_date = dcm.get('StudyDate', datetime.now().strftime('%d/%m/%Y'))
      if len(study_date)==8:
        study_date = f"{study_date[6:8]}/{study_date[4:6]}/{study_date[0:4]}"
      date_entry.insert(0, study_date)

      if (0x0010, 0x4000) in dcm:
        comment_entry.insert(0, dcm[0x00104000].value)
    except Exception as e:
      file_label.config(text=f"Błąd odczytu DICOM: {e}")
  else:
    im = loadImage(file_path)
    date_entry.insert(0, datetime.now().strftime('%d/%m/%Y'))

  if im is not None:
    showImage(im, 1)

#Stworzenie pola
def create_entry(parent, text):
  frame = ttk.Frame(parent)
  label = ttk.Label(frame, text=text)
  label.pack(side=tk.TOP, pady=(2, 0))
  entry = ttk.Entry(frame)
  entry.pack(side=tk.TOP, pady=(2, 0))
  frame.pack(padx=20)
  return entry


def image_fit(im, flag=False):
  img = Image.fromarray(rescale(im) if flag else im)
  img.thumbnail((300, 300))
  return ImageTk.PhotoImage(img)


def run_simulation():
  global im
  label_progress("")
  progress(0)
  if file_path is None:
    progress_bar_label.config(text="NIE WYBRANO PLIKU")
    return
  if not re.fullmatch(r"\d{2}/\d{2}/\d{4}", date_entry.get()):
    progress_bar_label.config(text="BŁĘDNA DATA")
    return
  if im is None:
    im = loadImage(file_path)

  #Obraz wejściowy
  showImage(im, 1)
  sinogram = radonTransform(im,
                               img_label=2,
                               numOfScans=int(360/int(step_entry.get())),
                               step=int(step_entry.get()),
                               emitterRange=int(emitter_range_entry.get()),
                               numOfDetectors=int(num_of_detectors_entry.get()))
  #Sinogram
  showImage(sinogram, 2)

  if filter_var.get() == 1:
    sinogram = filtr(sinogram, img_label=4)
    #Przefiltrowany sinogram
    showImage(sinogram, 4, flag=True)
  reconstructed = inverseRadonTransform(sinogram, im, img_label=5,
                                           numOfScans=int(360/int(step_entry.get())),
                                           step=int(step_entry.get()),
                                           emitterRange=int(emitter_range_entry.get()),
                                           numOfDetectors=int(num_of_detectors_entry.get()))
  #Obraz wyjściowy
  showImage(reconstructed, 5)
  if dicom_var.get() == 1:
    convertToDCM(reconstructed, name_entry.get(), id_entry.get(), date_entry.get(), comment_entry.get())
  else:
    Image.fromarray(reconstructed).convert("L").save("output.png")

#Wczytywanie obrazków
def showImage(im, img_l, flag=False):
  img_label_map = {1: input_photo_label, 2: sinogram_photo_label, 4: filtred_sinogram_photo_label, 5: output_photo_label}  # Mapowanie numerów na etykiety
  img_label = img_label_map.get(img_l)

  if img_label is None:
    print(f"Błąd: brak etykiety dla img_l = {img_l}")
    return

  img = image_fit(im, flag=flag)
  img_label.config(image=img)  # Aktualizacja obrazu
  img_label.image = img  # Zachowanie referencji do obrazu
  img_label.pack()
  root.update()


#Główny interfejs
root = ThemedTk(theme="ITFT1")
root.title("Projekt tomografu komputerowego")
root.geometry("1200x900")
file_path = None

#Kontener dla całej zawartości
main_container = ttk.Frame(root)
main_container.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)

#Górny panel kontrolny
control_panel = ttk.Frame(main_container)
control_panel.pack(side=tk.TOP, fill=tk.X, pady=(0, 10))

#Sekcja wyboru pliku
file_frame = ttk.Frame(control_panel)
file_frame.pack(side=tk.LEFT, padx=5)
choose_button = ttk.Button(file_frame, text="Wybierz plik", command=choose_file)
choose_button.pack()
file_label = ttk.Label(file_frame, anchor='center',  text="", width=50)
file_label.pack()

#Sekcja parametrów
params_frame = ttk.Frame(control_panel)
params_frame.pack(side=tk.LEFT, padx=20)

step_entry = ttk.Entry(params_frame, width=5)
step_entry.insert(0, '1')
step_entry.pack(side=tk.LEFT, padx=5)
ttk.Label(params_frame, text="Krok alfa").pack(side=tk.LEFT)

emitter_range_entry = ttk.Entry(params_frame, width=5)
emitter_range_entry.insert(0, '180')
emitter_range_entry.pack(side=tk.LEFT, padx=5)
ttk.Label(params_frame, text="Rozpiętość").pack(side=tk.LEFT)

num_of_detectors_entry = ttk.Entry(params_frame, width=5)
num_of_detectors_entry.insert(0, '180')
num_of_detectors_entry.pack(side=tk.LEFT, padx=5)
ttk.Label(params_frame, text="Detektory").pack(side=tk.LEFT)

#Sekcja opcji
options_frame = ttk.Frame(control_panel)
options_frame.pack(side=tk.LEFT, padx=20)

filter_var = tk.IntVar()
filter_cb = ttk.Checkbutton(options_frame, text="Filtrowanie", variable=filter_var)
filter_cb.pack(side=tk.LEFT, padx=5)

dicom_var = tk.IntVar()
dicom_cb = ttk.Checkbutton(options_frame, text="DICOM", variable=dicom_var, command=toggle_dicom_fields)
dicom_cb.pack(side=tk.LEFT, padx=5)

#Przycisk uruchomienia
run_button = ttk.Button(control_panel, text="Uruchom", command=run_simulation, width=20)
run_button.pack(side=tk.RIGHT)

#Sekcja danych DICOM
dicom_fields_frame = ttk.Frame(main_container)
dicom_fields_frame.pack(side=tk.TOP, fill=tk.X, pady=(5, 10))

#Pola DICOM obok siebie
name_frame = ttk.Frame(dicom_fields_frame)
name_frame.pack(side=tk.LEFT, padx=5)
ttk.Label(name_frame, text="Imię i nazwisko").pack()
name_entry = ttk.Entry(name_frame, width=25)
name_entry.pack()

id_frame = ttk.Frame(dicom_fields_frame)
id_frame.pack(side=tk.LEFT, padx=5)
ttk.Label(id_frame, text="ID pacjenta").pack()
id_entry = ttk.Entry(id_frame, width=15)
id_entry.pack()

date_frame = ttk.Frame(dicom_fields_frame)
date_frame.pack(side=tk.LEFT, padx=5)
ttk.Label(date_frame, text="Data badania").pack()
date_entry = ttk.Entry(date_frame, width=10)
date_entry.pack()

comment_frame = ttk.Frame(dicom_fields_frame)
comment_frame.pack(side=tk.LEFT, padx=5)
ttk.Label(comment_frame, text="Komentarz").pack()
comment_entry = ttk.Entry(comment_frame, width=30)
comment_entry.pack()

dicom_fields_frame.pack_forget()

#Pasek postępu
progress_frame = ttk.Frame(main_container)
progress_frame.pack(side=tk.TOP, fill=tk.X, pady=(0, 10))
progress_bar = ttk.Progressbar(progress_frame, orient='horizontal', mode='determinate', length=1000)
progress_bar.pack()
progress_bar_label = ttk.Label(progress_frame, text="", font=(30))
progress_bar_label.pack()

#Panel z obrazami
images_panel = ttk.Frame(main_container)
images_panel.pack(fill=tk.BOTH, expand=True)

#Kontener dla górnych obrazów
top_images_container = ttk.Frame(images_panel)
top_images_container.pack(side=tk.TOP, pady=(0, 20))
top_images = ttk.Frame(top_images_container)
top_images.pack()

#Górny rząd obrazów
input_photo = ttk.Frame(top_images)
input_photo.pack(side=tk.LEFT, padx=20)
input_photo_text = ttk.Label(input_photo, text="Wybrany obraz")
input_photo_text.pack()
input_photo_label = ttk.Label(input_photo)
input_photo_label.pack()

sinogram_photo = ttk.Frame(top_images)
sinogram_photo.pack(side=tk.LEFT, padx=20)
sinogram_photo_text = ttk.Label(sinogram_photo, text="Sinogram")
sinogram_photo_text.pack()
sinogram_photo_label = ttk.Label(sinogram_photo)
sinogram_photo_label.pack()

#Kontener dla dolnych obrazów
bottom_images_container = ttk.Frame(images_panel)
bottom_images_container.pack(side=tk.TOP)
bottom_images = ttk.Frame(bottom_images_container)
bottom_images.pack()

#Dolny rząd obrazów
filtred_sinogram_photo = ttk.Frame(bottom_images)
filtred_sinogram_photo.pack(side=tk.LEFT, padx=20)
filtred_sinogram_photo_text = ttk.Label(filtred_sinogram_photo, text="Przefiltrowany Sinogram")
filtred_sinogram_photo_text.pack()
filtred_sinogram_photo_label = ttk.Label(filtred_sinogram_photo)
filtred_sinogram_photo_label.pack()

output_photo = ttk.Frame(bottom_images)
output_photo.pack(side=tk.LEFT, padx=20)
output_photo_text = ttk.Label(output_photo, text="Obraz końcowy")
output_photo_text.pack()
output_photo_label = ttk.Label(output_photo)
output_photo_label.pack()

root.mainloop()