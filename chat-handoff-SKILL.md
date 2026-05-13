---
name: chat-handoff
description: >-
  Use when the user asks to summarize or hand off the current chat/context to a
  new session. Trigger phrases include: "суммаризируй чат", "суммаризируй
  диалог", "передай контекст", "сделай хандоф", "handoff", "выжимка диалога",
  "контекст для нового чата". Produce exactly one copy-paste markdown handoff
  block, with no trailing commentary.
---

# Chat Handoff — передача контекста в новую сессию

## Purpose

Turn the currently available conversation context into a clean, copy-pasteable markdown handoff for a new agent/session.

The handoff must help the next agent understand:

- what the user was trying to accomplish;
- what has already been decided, built, changed, rejected, or learned;
- what artifacts, files, links, commands, and constraints matter;
- what remains open;
- how to continue without asking the user to repeat known information.

This skill is for summarizing past context only. Do not insert a new task for the next agent inside the handoff.

## Trigger conditions

Use this skill when the user asks for any of the following or a close equivalent:

- summarize the current chat/dialogue;
- create a handoff for a new chat/session/agent;
- pass context forward;
- make a compact project summary;
- produce a copy-paste context block;
- phrases such as: `суммаризируй чат`, `передай контекст`, `сделай хандоф`, `handoff`, `выжимка диалога`, `контекст для нового чата`.

Do not require the chat to be technical. This skill works for research, coding, writing, planning, design, decision-making, and mixed conversations.

## Non-goals

Do not use this skill to:

- answer a new substantive question;
- continue implementation work;
- rewrite the user's project unless the user explicitly asks for that separately;
- invent missing context;
- produce a chronological transcript;
- create a marketing-style recap;
- summarize external material that has not been provided or discussed.

## Core rules

1. **Summarize only available context.** Use the conversation, visible tool outputs, uploaded files, created files, and artifacts that are available in the current session. If something is not available, say so only if it matters.
2. **Do not hallucinate.** Every fact, file path, decision, reason, quote, and preference must be grounded in the chat.
3. **No new task inside the handoff.** The handoff is read-only context. The user will write the next instruction separately.
4. **No empty sections.** Include only sections that contain real information.
5. **Concrete beats generic.** Prefer `src/db.ts stores SQLite setup` over `we discussed code`.
6. **Decisions need reasons only when reasons were actually stated.** If the chat contains no reason, list the decision without inventing a `because`.
7. **Preserve rejected alternatives.** If options were considered and rejected, include them so the next agent does not repeat them.
8. **User preferences require evidence.** Only mention style/work preferences when the user explicitly stated or demonstrated them in the chat.
9. **Language follows the recent conversation.** Write the handoff in the language used by the user in the last few substantive messages. Keep technical terms in their original language when useful.
10. **One final answer only.** Output one short instruction sentence, then one fenced markdown block, and nothing after the block.

## Workflow

### 1. Collect context silently

Review the available session context and identify:

- initial goal: what the user wanted at the start;
- current goal: how the goal changed or became more precise;
- key facts, findings, and conclusions;
- decisions and the stated reasons for them;
- artifacts: files created/edited, paths, generated docs, commands, code snippets, links;
- problems encountered and how they were resolved;
- external sources, references, datasets, screenshots, or documents used;
- alternatives considered and rejected;
- terminology introduced in the chat;
- explicit constraints, assumptions, requirements, and non-goals;
- user style preferences that are evidenced by the conversation;
- open questions and likely next steps.

Do this silently. Do not narrate the collection process.

### 2. Detect the conversation type

Choose the dominant type internally. Do not announce it.

| Type | Signs | Output implication |
|---|---|---|
| `research` | reading sources, comparing ideas, extracting conclusions | emphasize findings, sources, open questions |
| `project/build` | code, files, commands, implementation, debugging | emphasize artifacts, paths, decisions, known issues |
| `creative` | text, design, naming, content iterations | emphasize direction, chosen style, rejected versions |
| `decision` | choosing between options or strategies | emphasize criteria, alternatives, decisions, tradeoffs |
| `mixed` | multiple unrelated tracks | split by topic so each topic is self-contained |

### 3. Ask clarification only if truly blocking

Usually skip questions and produce the handoff immediately.

Ask at most 1–3 concise questions only when a real ambiguity would make the handoff misleading, for example:

- which of several unrelated contexts should be handed off;
- whether to include confidential/sensitive details;
- whether the user wants a compact or expanded handoff and the chat is very long.

If the ambiguity is minor, make the safest reasonable choice and continue.

### 4. Build the handoff

Output format:

1. First line outside the block:
   `Скопируй блок ниже и вставь первым сообщением в новый чат.`
2. Then one fenced markdown block.
3. Nothing after the block.

Use a four-backtick outer fence when the handoff itself may contain triple-backtick code snippets.

## Output templates

### Compact template

Use for short chats, low complexity, or when there are fewer than about 10 substantive messages.

````markdown
# Контекст из прошлого чата

> Это передача контекста из предыдущей сессии. Прочитай, подтверди коротким «принял», и жди следующего сообщения пользователя.

**TL;DR:** <one sentence with the main result or state>

## О чём был чат
<1–2 sentences>

## Главное
- <key fact, decision, or result>
- <key fact, decision, or result>

## Открытые вопросы / следующий шаг
- <what remains unresolved or the next natural continuation>
````

### Standard template

Use for most chats.

````markdown
# Контекст из прошлого чата

> Это передача контекста из предыдущей сессии. Прочитай, подтверди коротким «принял», и жди следующего сообщения пользователя.

**TL;DR:** <one sentence that makes the whole handoff understandable without reading further>

## О чём был чат
<1–3 sentences explaining the topic and purpose>

## Главное, что выяснили / сделали
- <specific finding, result, or completed work>
- <specific finding, result, or completed work>

## Ключевые решения и почему
- **<decision>** — <reason stated in the chat, or omit the reason if none was stated>
- **<decision>** — <reason stated in the chat, or omit the reason if none was stated>

## Источники и материалы
- <file/link/document/source> — <why it matters>

## Стиль работы и предпочтения пользователя
- <preference with clear evidence from the chat>

## Открытые вопросы / куда двигаться дальше
- <unresolved issue or next continuation point>
````

### Mixed-topic template

Use when the chat contains two or more unrelated tracks.

````markdown
# Контекст из прошлого чата

> Это передача контекста из предыдущей сессии. Прочитай, подтверди коротким «принял», и жди следующего сообщения пользователя.

**TL;DR:** Чат был о <N> темах: <topic 1>, <topic 2>, <topic 3>. Ниже каждая тема вынесена отдельно.

## О чём был чат
<Briefly explain why this is a mixed handoff and what the independent topics are.>

## Тема 1: <name>
- **Главное:** <facts/results>
- **Решения:** <decisions and stated reasons>
- **Артефакты / источники:** <files, links, materials, if any>
- **Открытые вопросы:** <what remains unresolved, if any>

## Тема 2: <name>
- **Главное:** <facts/results>
- **Решения:** <decisions and stated reasons>
- **Артефакты / источники:** <files, links, materials, if any>
- **Открытые вопросы:** <what remains unresolved, if any>

## Стиль работы и предпочтения пользователя
- <cross-topic preferences with evidence from the chat>
````

## Optional sections

Add these only when there is real content. Place them before `Открытые вопросы / куда двигаться дальше` in a logical order.

| Section heading | Add when |
|---|---|
| `## Артефакты (файлы и пути)` | files, docs, code, generated assets, or paths were created/changed/discussed |
| `## Команды и проверки` | shell commands, tests, build steps, migrations, or verification commands matter |
| `## Решённые проблемы` | there was debugging, an obstacle, or a misunderstanding that was resolved |
| `## Рассмотренные альтернативы` | options were compared, rejected, or deferred |
| `## Терминология` | custom terms, names, abbreviations, or domain-specific vocabulary were introduced |
| `## Ограничения и допущения` | explicit constraints, assumptions, requirements, or non-goals were stated |
| `## Ключевые цитаты` | exact wording from a source or user instruction is important |
| `## Риски / что не сломать` | future work could accidentally undo important decisions or constraints |
| `## Состояние репозитория` | repo state, branch, changed files, test state, or uncommitted work matters |

## Style guidance

- Write densely but clearly.
- Prefer bullets over long paragraphs.
- Keep the TL;DR to one sentence.
- Do not include greetings, apologies, or meta-commentary.
- Do not use emoji in section headings.
- Do not say “we first did X, then Y” unless sequence is essential.
- For technical chats, preserve exact filenames, commands, package names, errors, and versions when available.
- For research chats, distinguish facts, interpretations, and open questions.
- For decision chats, preserve criteria and rejected alternatives.
- For creative chats, preserve the chosen direction, tone, examples, and rejected styles.

## Handling missing or uncertain context

If the current agent cannot see earlier parts of the conversation, do not pretend otherwise.

Use a concise note inside the handoff only if necessary:

- `Контекст до этого сообщения недоступен; выжимка основана только на видимой части диалога.`
- `Файл упоминался в чате, но его содержимое недоступно в текущей сессии.`

Do not overuse uncertainty notes. Include them only when they prevent a misleading handoff.

## Final quality checklist

Before sending, verify silently:

1. The handoff contains only past context, not a new instruction.
2. The TL;DR is useful on its own.
3. Every included section has real content.
4. All facts are grounded in visible context.
5. File paths, commands, URLs, names, and versions are exact when known.
6. Decisions include reasons only when reasons were stated.
7. Rejected alternatives are included when they matter.
8. User preferences are included only when evidenced.
9. The language matches the user's recent language.
10. The response has one instruction sentence, one fenced markdown block, and no text after it.
