import numpy as np
import streamlit as st
import os

# punc = '''!()-[]{};:'"\,<>./?@#$%^&*_~'''
characters = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'

def remove_punc(l):
    split_l = l.split()
    new_l = []
    for word in split_l:
        no_punc = ''
        for char in word:
            if char.upper() in characters:
                no_punc += char
        if no_punc:
            new_l.append(no_punc.lower())
    return new_l

def poem_to_words(txt_file):
    saved_words = []
    with open(txt_file, 'r') as f:
        lines = f.readlines()
        for l in lines:
            # check that line has data
            if len(l.strip()) > 0:
                # split and remove punctuation
                split_l = remove_punc(l)
                for word in split_l:
                    saved_words.append(word)
    return saved_words

def poem_to_list(txt_file, user_file):
    saved_lines = []

    if not user_file:
        with open(txt_file, 'r') as f:
            lines = f.readlines()
            for l in lines:
                # check that line has data
                if len(l.strip()) > 0:
                    # split and remove punctuation
                    split_l = remove_punc(l)
                    saved_lines.append(split_l)
    
    else:
        for l in txt_file.split('\n'):
            # check that line has data
            if len(l.strip()) > 0:
                # split and remove punctuation
                split_l = remove_punc(l)
                saved_lines.append(split_l)
    
    return saved_lines
               

def generate_models(poet_list):

    # don't judge. i realized late that i needed int2word and word2int, so this was quicker than thinking
    word2int = {}
    idx = 0
    for line in poet_list:
        for word in line:
            if word not in word2int:
                word2int[word] = idx
                idx += 1

    pi, A1, A2 = {}, {}, {}

    for line in poet_list:
        for idx, word in enumerate(line):
            if idx == 0:
                # store word in pi
                if word not in pi:
                    pi[word2int[word]] = 0
                pi[word2int[word]] += 1.
            if idx < len(line) - 1:
                # store word in A1
                if word2int[word] not in A1:
                    A1[word2int[word]] = {}
                word_plus1 = line[idx+1]
                if word_plus1 not in A1[word2int[word]]:
                    A1[word2int[word]][word2int[word_plus1]] = 0
                A1[word2int[word]][word2int[word_plus1]] += 1.
            if idx < len(line) - 2:
                # store word in A2:
                if word2int[word] not in A2:
                    A2[word2int[word]] = {}
                if word2int[word_plus1] not in A2[word2int[word]]:
                    A2[word2int[word]][word2int[word_plus1]] = {}
                word_plus2 = line[idx+2]
                if word2int[word_plus2] not in A2[word2int[word]][word2int[word_plus1]]:
                    A2[word2int[word]][word2int[word_plus1]][word2int[word_plus2]] = 0
                A2[word2int[word]][word2int[word_plus1]][word2int[word_plus2]] += 1.
    
    return pi, A1, A2, word2int

def norm_dict_values(d):
    if len(d) == 1:
        for k in d.keys():
            return {k: 1.}
    cum_sum = sum(d.values())
    probs = {}
    for k,v in d.items():
        probs[k] = v / cum_sum
    
    return probs

def norm_models(pi, A1, A2):
    A2_norm = {}
    A1_norm = {}

    for k1 in A2:
        if k1 not in A2_norm:
            A2_norm[k1] = {}
        for k2 in A2[k1]:
            A2_norm[k1][k2] = norm_dict_values(A2[k1][k2])

    for k1 in A1:
        A1_norm[k1] = norm_dict_values(A1[k1])

    pi_norm = norm_dict_values(pi)

    return pi_norm, A1_norm, A2_norm

def models_from_poet(poet):
    pi, A1, A2, word2int = generate_models(poet)
    pi_norm, A1_norm, A2_norm = norm_models(pi, A1, A2)

    return pi_norm, A1_norm, A2_norm, word2int


# time to generate
# get first word from pi_norm dist
# get second word from A1_norm[word - 1] dist
# get rest of line from A2_norm[word - 2][word - 1] dist

# test with a 10 word sentence

def generate_sentence(max_len, pi_norm, A1_norm, A2_norm, word2int):

    int2word = {v: k for k, v in word2int.items()}
    
    for idx in range(max_len):
        if idx == 0:
            first_word = np.random.choice(list(pi_norm.keys()), p=list(pi_norm.values()))
            int_sentence = [first_word]
            word_sentence = int2word[first_word]
        elif idx == 1:
            if first_word not in A1_norm:
                return f'{word_sentence}.'
            second_word = np.random.choice(list(A1_norm[first_word].keys()), p=list(A1_norm[first_word].values()))
            int_sentence.append(second_word)
            word_sentence += f' {int2word[second_word]}'
        else:
            word_minus2 = int_sentence[idx-2]
            word_minus1 = int_sentence[idx-1]

            if word_minus2 not in A2_norm:
                return f'{word_sentence}.'
                # word_sentence += '.'
            if word_minus1 not in A2_norm[word_minus2]:
                return f'{word_sentence}.'
                    
            else:
                choices = list(A2_norm[word_minus2][word_minus1].keys())
                probs = list(A2_norm[word_minus2][word_minus1].values())
                following_word = np.random.choice(choices, p=probs)
                
                int_sentence.append(following_word)
                word_sentence += f' {int2word[following_word]}'

    return word_sentence
    

def generate_poem(models, num_lines, max_line_len):
    sentences = []
    for _ in range(num_lines):
        sentence = generate_sentence(max_line_len, *models)
        # st.markdown(sentence.capitalize())
        sentences.append(sentence.capitalize())
    st.session_state.poem = sentences

def reset_state():
    st.session_state.poem = None

st.set_page_config(
    page_title='Poetry Generator',
    page_icon=':scroll:',
    layout='wide'
)

st.title(' :scroll: Poetry Generator :scroll: ')
st.markdown('**Set options in the sidebar**')

if 'poem' not in st.session_state:
    st.session_state.poem = None

# sidebar
st.sidebar.header('Select options:')
poet_selector = st.sidebar.selectbox('Select poet:', ('Shel Silverstein', 'Edgar Allan Poe', 'Robert Frost', 'Learn from my poems'), on_change=reset_state)
num_lines = st.sidebar.slider('Number of lines', min_value=1, max_value=8, value=4, step=1)
max_line_len = st.sidebar.slider('Maximum number of words per sentence', min_value=2, max_value=14, value=8, step=1)

if st.session_state.poem:
    for line in st.session_state.poem:
        st.markdown(line)

if poet_selector:
    if poet_selector != 'Learn from my poems':
        author2file = {
            'Shel Silverstein': 'shel_silverstein',
            'Robert Frost': 'robert_frost',
            'Edgar Allan Poe': 'edgar_allan_poe'
        }

        poet_list = poem_to_list(f'./poem_generator/poem_inputs/{author2file[poet_selector]}.txt')
        # poet_list = poem_to_list(f'./poem_inputs/{author2file[poet_selector]}.txt', user_file=False)
        models = models_from_poet(poet_list)
        
    else: 
        user_poems = st.text_area('Upload poems to learn from')
        # st.markdown(user_poems)
        poet_list = poem_to_list(user_poems, user_file=True)
        models = models_from_poet(poet_list)
    
    st.button('Generate poem!', on_click=generate_poem, args=[models, num_lines, max_line_len])


