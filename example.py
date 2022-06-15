from reedsolomon import ReedSolomonCodec

codec = ReedSolomonCodec(9) #10 parnostnih bitov

message = "Zivjo Damjan, kako si?" #sporočilo je lahko seznam vrednosti 0-255 ali pa niz ASCII simbolov (brez šumnikov)
encoded = codec.encode(message) #obstaja tudi parameter 'return_str', ki potem povzroči da encode metoda vrne niz
errors = codec.erasure_sim(encoded, [1, 6, 3, 10, 7]) #metoda simulira izbris črk/vrednosti na teh položajih in vrne sporočilo z napakami

decoded = codec.decode_erasures(errors, [1, 6, 3, 10, 7], return_str=True) #pomembno je, da so indeksi izbrisov točni, drugače sporočilo ne bo enako začetnemu
print(decoded) #vrnjeno je bilo zakodirano sporočilo z dodatnimi 9-imi znaki